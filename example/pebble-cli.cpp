/*!
  \file pebble-cli.cpp
  \brief Classify a pebble stored in an STL file.
*/

#include <boost/graph/properties.hpp>
//#define PRINT 1

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <forward_list>
#include <cstring>
#include <boost/container/static_vector.hpp>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <boost/math/constants/constants.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>

#include <pebble/utility.hpp>
#include <pebble/convex_hull.hpp>
#include <pebble/face_center.hpp>
#include <pebble/edge_type.hpp>
#include <pebble/edge_saddle.hpp>
#include <pebble/vertex_dividing.hpp>
#include <pebble/contour_line.hpp>
#include <pebble/ascending_curve.hpp>
#include <pebble/descending_curve.hpp>
#include <pebble/intersections.hpp>
#include <pebble/sort_saddles.hpp>
#include <pebble/cancel_saddle.hpp>
#include <pebble/create_reeb.hpp>
#include <pebble/encode_graph.hpp>
#include <pebble/create_morse_smale.hpp>
#include <pebble/get_saddles.hpp>
#include <pebble/create_master_graph.hpp>
#include <pebble/distance_from_centroid.hpp>
#include <pebble/axes.hpp>
#include <pebble/quasi_dual.hpp>
#include <pebble/measures.hpp>
#include <pebble/fit_paraboloid.hpp>
#include <pebble/master_weight.hpp>

#include "types.hpp"

using CGAL::Polygon_mesh_processing::compute_face_normals;

/*!
  \brief Print the help message on the given output stream.
  \param stream the help text will be outputted to this stream
 */
void print_help(std::ostream &stream) {
  stream << "Usage: pebble-cli FILE STABLE UNSTABLE\n"
	 << "Calculate shape descriptors of the pebble in FILE and in the primary equilibrium class {STABLE,UNSTABLE}.\n\n"
	 << "FILE must be in STereoLitography (STL) format.\n\n"
	 << "The following values are printed in a semicolon separated list:\n"
	 << "  - filename\n"
	 << "  - STABLE and UNSTABLE parameters\n"
	 << "  - axes of the approximating ellipsoid (a>b>c) and the axis ratios c/b and b/a\n"
	 << "  - Zingg class of the pebble (disc, sphere, blade or rod)\n"
	 << "  - volume (V) and surface area (A) of the mesh\n"
	 << "  - isoperimetric ratio of the mesh (36*pi*V^2/A^3)\n"
	 << "  - computed primary equilibrium class, approximating {STABLE,UNSTABLE} the best\n"
	 << "  - an alphanumerical encoding of the master, Reeb and Morse-Smale graphs, and also the quasi dual of the Morse-Smale graph (the encodings of two graphs are equal if and only if they are isomorphic)\n";
}

template <typename EdgeTypeMap>
struct AscendingEdge {
  EdgeTypeMap m_type;
  AscendingEdge() {}
  AscendingEdge(EdgeTypeMap edge_type) : m_type(edge_type) {}
  template <typename Edge>
  bool operator()(const Edge &edge) const {
    return get(m_type, edge) == pebble::EdgeType::ASCENDING;
  }
};

template <typename VertexTypeMap>
struct TypedVertex {
  VertexTypeMap m_poly_type;
  bool m_use_poly;
  VertexTypeMap m_para_type;
  pebble::VertexType m_target;
  TypedVertex() {}
  TypedVertex(VertexTypeMap poly_type,
	      bool use_poly,
	      VertexTypeMap para_type,
	      pebble::VertexType target)
    : m_poly_type(poly_type),
      m_use_poly(use_poly),
      m_para_type(para_type),
      m_target(target) {}
  template <typename Vertex>
  bool operator()(const Vertex &vertex) const {
    return get(m_para_type, vertex) == m_target ||
      (!m_use_poly || get(m_poly_type, vertex) == pebble::VertexType::SADDLE) &&
      get(m_para_type, vertex) == pebble::VertexType::INTERSECTION;
  }
};

using EdgeTypeMap =
    boost::property_map<UndirectedMaster, pebble::EdgeType EdgeData::*>::type;
using VertexTypeMap =
    boost::property_map<UndirectedMaster, pebble::VertexType VertexData::*>::type;
using FilteredGraph = boost::filtered_graph<UndirectedMaster, AscendingEdge<EdgeTypeMap>,
                                            TypedVertex<VertexTypeMap>>;

int main(int argc, char **argv) {
  using namespace boost::math::constants;
  // parse options
  if (argc == 2 && std::strncmp(argv[1], "-h", 3) == 0) {
    print_help(std::cout);
    return EXIT_SUCCESS;
  }
  if (argc < 4) {
    print_help(std::cerr);
    return EXIT_FAILURE;
  }
  Mesh mesh;
  std::unordered_map<FaceHandle, bool> original_storage;
  auto face_original_pm = boost::make_assoc_property_map(original_storage);
  const bool read = pebble::convex_hull(argv[1], mesh, face_original_pm);
  if (!read || mesh.num_vertices() == 0) {
    std::cerr << "pebble-cli: could not read file\n"
	      << "Try 'pebble-cli -h' for more information.\n";
    return EXIT_FAILURE;
  }

  std::cerr << std::fixed;
#ifdef PRINT
  std::cerr << "color cyan\n";
  for (const FaceHandle &face : mesh.faces()) {
    const HalfedgeHandle he = halfedge(face, mesh);
    if (face_original_pm[face])
      std::cerr << "triangle ("
		<< mesh.point(mesh.source(he)) << ") ("
		<< mesh.point(mesh.target(he)) << ") ("
		<< mesh.point(mesh.target(mesh.next(he))) << ")\n";
  }
  std::cerr << "color yellow\n";
  for (const FaceHandle &face : mesh.faces()) {
    const HalfedgeHandle he = halfedge(face, mesh);
    if (!face_original_pm[face])
      std::cerr << "triangle ("
		<< mesh.point(mesh.source(he)) << ") ("
		<< mesh.point(mesh.target(he)) << ") ("
		<< mesh.point(mesh.target(mesh.next(he))) << ")\n";
  }
#endif

  const int stable_target = atoi(argv[2]);
  const int unstable_target = atoi(argv[3]);
  if (stable_target < 1 || unstable_target < 1) {
    std::cerr << "pebble-cli: could not parse primary equilibrium class\n"
	      << "Try 'pebble-cli -h' for more information.\n";
    return EXIT_FAILURE;
  }

  const Point centroid = CGAL::centroid(mesh.points().begin(),
					mesh.points().end());

  // prepare some commonly used properties
  const auto face_normals_pm = mesh.add_property_map<FaceHandle, Vector>("f:normal").first;
  compute_face_normals(mesh, face_normals_pm);

  const auto face_center_pm =
    mesh.add_property_map<FaceHandle, Point>("f:center").first;
  const auto face_barycentric_pm =
    mesh.add_property_map<FaceHandle, Barycentric>("f:barycentric").first;
  const auto is_followed_pm =
    mesh.add_property_map<EdgeHandle, bool>("e:followed", false).first;
  const auto dividing =
    mesh.add_property_map<HalfedgeHandle, bool>("h:dividing", false).first;
  const auto edge_minimum_pm =
    mesh.add_property_map<EdgeHandle, Point>("e:minimum").first;
  const auto edge_minimum_barycentric =
    mesh.add_property_map<EdgeHandle, std::array<Scalar, 2>>("e:minimum_barycentric").first;

  // face and edge minimum
  pebble::face_center(mesh, centroid, face_center_pm, face_barycentric_pm, mesh.points());
  pebble::edge_saddle(mesh, centroid, edge_minimum_pm, edge_minimum_barycentric, mesh.points());

  // classify edges and vertices based on the position of the face minimum
  pebble::edge_type(mesh, face_barycentric_pm, is_followed_pm);
  pebble::vertex_dividing(mesh, face_barycentric_pm, dividing);

  // find all the saddle points
  std::vector<EdgeHandle> edges_with_saddle;
  pebble::get_saddles(mesh,
                      is_followed_pm,
                      edge_minimum_barycentric,
                      std::back_inserter(edges_with_saddle));

  // trace contour lines
  auto halfedge_contours_pm =
      mesh.add_property_map<HalfedgeHandle,
                            std::vector<pebble::HalfedgeContour<HalfedgeHandle,
                                                                Point>>>("h:contour").first;
  auto face_contours_pm =
      mesh.add_property_map<FaceHandle,
                            std::vector<pebble::FaceContour<HalfedgeHandle,
                                                            Point, Scalar>>>("f:contour").first;
  pebble::contour_lines(mesh,
                        centroid,
			mesh.points(),
                        edge_minimum_pm,
                        edge_minimum_barycentric,
                        face_center_pm,
                        face_normals_pm,
                        edges_with_saddle.cbegin(),
                        edges_with_saddle.cend(),
                        halfedge_contours_pm,
                        face_contours_pm);

  //trace ascending curves
  auto ascending_segments_map =
      mesh.add_property_map<HalfedgeHandle,
                            std::vector<pebble::ascending_segment<Mesh,
                                                                  Segment>>>("h:ascending_segments").first;
  auto unstable_map =
    mesh.add_property_map<HalfedgeHandle, VertexHandle>("h:unstable").first;
  pebble::ascending_curves(mesh,
                           centroid,
			   mesh.points(),
                           face_center_pm,
                           edge_minimum_pm,
                           edge_minimum_barycentric,
                           is_followed_pm,
                           dividing,
                           edges_with_saddle.cbegin(),
                           edges_with_saddle.cend(),
                           ascending_segments_map,
                           unstable_map);

  // trace descending curves
  auto descending_segments_map =
      mesh.add_property_map<HalfedgeHandle,
                            std::forward_list<pebble::ascending_segment<Mesh,
                                                                        Segment>>>("h:descending_segments").first;
  auto stable_map =
    mesh.add_property_map<HalfedgeHandle, FaceHandle>("h:stable").first;
  pebble::descending_curves(mesh,
                            centroid,
			    mesh.points(),
                            edge_minimum_pm,
                            face_center_pm,
                            face_barycentric_pm,
                            edges_with_saddle.cbegin(),
                            edges_with_saddle.cend(),
                            descending_segments_map,
                            stable_map);

  // find intersection points
  auto ascending_intersections_map =
      mesh.add_property_map<HalfedgeHandle,
                            std::vector<pebble::Intersection<HalfedgeHandle,
                                                             Point>>>("h:ascending_intersections").first;
  pebble::intersections(mesh,
                        centroid,
                        face_center_pm,
                        face_normals_pm,
                        face_contours_pm,
                        halfedge_contours_pm,
                        ascending_segments_map,
                        edges_with_saddle.cbegin(),
                        edges_with_saddle.cend(),
                        ascending_intersections_map);
  auto descending_intersections_map =
    mesh.add_property_map<HalfedgeHandle,
                          std::vector<pebble::Intersection<HalfedgeHandle,
                                                           Point>>>("h:descending_intersections").first;
  pebble::intersections(mesh,
                        centroid,
                        face_center_pm,
                        face_normals_pm,
                        face_contours_pm,
                        halfedge_contours_pm,
                        descending_segments_map,
                        edges_with_saddle.cbegin(),
                        edges_with_saddle.cend(),
                        descending_intersections_map);

  // build the master graph
  Master master;
  pebble::create_master_graph(mesh,
                              mesh.points(),
                              edge_minimum_pm,
                              face_center_pm,
                              edges_with_saddle.cbegin(),
                              edges_with_saddle.cend(),
                              ascending_intersections_map,
                              unstable_map,
                              descending_intersections_map,
                              stable_map,
                              master,
                              get(&VertexData::point, master),
                              get(&VertexData::type, master),
                              get(&VertexData::simplex, master),
                              get(&EdgeData::type, master));

  unsigned int index = 0;
  for (const MasterVertex vertex : boost::make_iterator_range(vertices(master)))
    master[vertex].id = index++;

  pebble::fit_paraboloid(mesh,
			 mesh.points(),
			 face_original_pm,
			 centroid,
			 master,
			 get(&VertexData::point, master),
			 get(&VertexData::type, master),
			 get(&VertexData::simplex, master),
			 get(&VertexData::para_type, master));

  int stable_target2 = 0;
  int unstable_target2 = 0;
  for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master))) {
    switch (master[vertex].para_type) {
    case pebble::VertexType::MIN:
      stable_target2++;
#ifdef PRINT
      std::cerr << "color green\npoint ("
		<< master[vertex].point << ")\n";
#endif
      break;
    case pebble::VertexType::MAX:
      unstable_target2++;
#ifdef PRINT
      std::cerr << "color red\npoint ("
		<< master[vertex].point << ")\n";
#endif
      break;
    }
  }

  // UndirectedMaster undir;
  // boost::copy_graph(master, undir, boost::vertex_index_map(get(&VertexData::id, master)));

  // AscendingEdge<EdgeTypeMap> edge_filter(get(&EdgeData::type, undir));
  // TypedVertex<VertexTypeMap> stable_filter(get(&VertexData::type, undir),
  // 					   false,
  // 					   get(&VertexData::para_type, undir),
  // 					   pebble::VertexType::MIN);
  // TypedVertex<VertexTypeMap> unstable_filter(get(&VertexData::type, undir),
  // 					     false,
  // 					     get(&VertexData::para_type, undir),
  // 					     pebble::VertexType::MAX);
  // FilteredGraph master_stable(undir, edge_filter, stable_filter);
  // FilteredGraph master_unstable(undir, edge_filter, unstable_filter);

  // for (const auto &vertex : boost::make_iterator_range(vertices(master_stable)))
  //   master_stable[vertex].color = boost::default_color_type();
  // boost::connected_components(master_stable,
  // 			      get(&VertexData::flock, master_stable),
  // 			      boost::color_map(get(&VertexData::color,
  // 						   master_stable)));
  // std::set<unsigned int> flocks_with_stable;
  // for (const auto &vertex : boost::make_iterator_range(vertices(master_stable)))
  //   if (master_stable[vertex].para_type == pebble::VertexType::MIN)
  //     flocks_with_stable.insert(master_stable[vertex].flock);
  // const int stable_target2 = flocks_with_stable.size();

  // for (const auto &vertex : boost::make_iterator_range(vertices(master_unstable)))
  //   master_unstable[vertex].color = boost::default_color_type();
  // boost::connected_components(master_unstable,
  // 			      get(&VertexData::flock, master_unstable),
  // 			      boost::color_map(get(&VertexData::color,
  // 						   master_unstable)));
  // std::set<unsigned int> flocks_with_unstable;
  // for (const auto &vertex : boost::make_iterator_range(vertices(master_unstable)))
  //   if (master_unstable[vertex].para_type == pebble::VertexType::MAX)
  //     flocks_with_unstable.insert(master_unstable[vertex].flock);
  // const int unstable_target2 = flocks_with_unstable.size();

  // set the hierarchy of saddles in the master graph
  pebble::distance_from_centroid(master,
                                 centroid,
                                 get(&VertexData::point, master),
                                 get(&VertexData::distance, master));
  std::vector<MasterVertex> saddles_sorted;
  pebble::sort_saddles(master,
                       get(&VertexData::type, master),
                       get(&VertexData::distance, master),
                       get(&EdgeData::type, master),
                       std::back_inserter(saddles_sorted));

  // {S, U} pair before canceling the saddle at the same index,
  // the last saddle will not be cancelled
  int stable = pebble::count_type(master,
				  get(&VertexData::type, master),
				  pebble::VertexType::MIN);
  int unstable = pebble::count_type(master,
				    get(&VertexData::type, master),
				    pebble::VertexType::MAX);
  std::vector<std::pair<int, int>> su(1, std::make_pair(stable, unstable));
  std::vector<std::string> masters, reebs, morse_smales, duals;

  // walk the hierarchy of master graphs
  // and generate the Reeb and Morse-Smale graphs
  {
    masters.push_back(pebble::encode_graph(master));
    {
      Master reeb = pebble::create_reeb(master, get(&EdgeData::type, master));
      unsigned int index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(reeb)))
	reeb[vertex].id = index++;
      reebs.push_back(pebble::encode_graph(reeb));
    }
    {
      Master ms = pebble::create_morse_smale(master,
                                             get(&VertexData::type, master),
                                             get(&EdgeData::type, master));
      unsigned int index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(ms)))
	ms[vertex].id = index++;
      morse_smales.push_back(pebble::encode_graph(ms));
      
      Master dual = pebble::quasi_dual(ms, get(&VertexData::type, master));
      index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(dual)))
	dual[vertex].id = index++;
      duals.push_back(pebble::encode_graph(dual));
    }
  }
  for (const MasterVertex &saddle : saddles_sorted) {
    pebble::cancel_saddle(master,
                          get(&VertexData::type, master),
                          get(&VertexData::distance, master),
                          get(&EdgeData::type, master),
                          saddle);

    unsigned int index = 0;
    for (const MasterVertex vertex : boost::make_iterator_range(vertices(master)))
      master[vertex].id = index++;

    stable = pebble::count_type(master,
				get(&VertexData::type, master),
				pebble::VertexType::MIN);
    unstable = pebble::count_type(master,
				  get(&VertexData::type, master),
				  pebble::VertexType::MAX);

    su.emplace_back(stable, unstable);
    masters.push_back(pebble::encode_graph(master));
    {
      Master reeb = pebble::create_reeb(master, get(&EdgeData::type, master));
      unsigned int index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(reeb)))
	reeb[vertex].id = index++;
      reebs.push_back(pebble::encode_graph(reeb));
    }
    {
      Master ms = pebble::create_morse_smale(master,
                                             get(&VertexData::type, master),
                                             get(&EdgeData::type, master));
      unsigned int index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(ms)))
	ms[vertex].id = index++;
      morse_smales.push_back(pebble::encode_graph(ms));

      Master dual = pebble::quasi_dual(ms, get(&VertexData::type, master));
      index = 0;
      for (const MasterVertex vertex : boost::make_iterator_range(vertices(dual)))
	dual[vertex].id = index++;
      duals.push_back(pebble::encode_graph(dual));
    }
  }

  // find the {S, U} pair closest to the target
  std::vector<int> difference;
  difference.reserve(su.size());
  for (const auto [s,u] : su)
    difference.push_back(abs(u - unstable_target) + abs(s - stable_target));
  const std::size_t min_index = std::distance(difference.begin(),
					      min_element(difference.begin(),
							  difference.end()));
  std::vector<int> difference2;
  difference2.reserve(su.size());
  for (const auto [s,u] : su)
    difference2.push_back(abs(u - unstable_target2));
  const std::size_t min_index2 = std::distance(difference2.begin(),
					       min_element(difference2.begin(),
							   difference2.end()));

  const std::array<Vector, 3> abc = pebble::axes(mesh, mesh.points());
  const Scalar a = sqrt(abc[0].squared_length());
  const Scalar b = sqrt(abc[1].squared_length());
  const Scalar c = sqrt(abc[2].squared_length());
  const Scalar y1 = c / b;
  const Scalar y2 = b / a;
  const Scalar vol = pebble::volume(mesh, mesh.points(), centroid);
  const Scalar area = pebble::surface_area(mesh, mesh.points());

  std::cout << std::filesystem::path(argv[1]).stem().string() << ';'
	    << stable_target << ';'
	    << unstable_target << ';'
	    << stable_target2 << ';'
	    << unstable_target2 << ';'
	    << a << ';'
	    << b << ';'
	    << c << ';'
	    << y1 << ';'
	    << y2 << ';'
	    << ((y1 < 2.0 / 3.0)
		? ((y2 > 2.0 / 3.0) ? "disc" : "blade")
		: ((y2 > 2.0 / 3.0) ? "sphere" : "rod")) << ';'
	    << vol << ';'
	    << area << ';'
	    << (36 * pi<Scalar>() * pow(vol, 2.0) / pow(area, 3.0)) << ';'
	    << su[min_index].first << ';'
	    << su[min_index].second << ';'
	    << su[min_index2].first << ';'
	    << su[min_index2].second << ';'
	    << masters[min_index] << ';'
	    << reebs[min_index] << ';'
	    << morse_smales[min_index] << ';'
	    << duals[min_index];

  return EXIT_SUCCESS;
}
