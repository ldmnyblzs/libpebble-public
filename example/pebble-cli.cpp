/*!
  \file pebble-cli.cpp
  \brief Classify a pebble stored in an STL file.
*/

#include <filesystem>
#include <fstream>
#include <forward_list>
#include <cstring>
#include <boost/container/static_vector.hpp>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <boost/math/constants/constants.hpp>

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

#include "types.hpp"

using CGAL::Polygon_mesh_processing::compute_face_normals;

int count_type(const Master &master,
	       const pebble::VertexType &type) {
  // TODO: this is not very optimal
  int result = 0;
  for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
    if (master[vertex].type == type)
      ++result;
  return result;
}

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
  const bool read = pebble::convex_hull(argv[1], mesh);
  if (!read || mesh.num_vertices() == 0) {
    std::cerr << "pebble-cli: could not read file\n"
	      << "Try 'pebble-cli -h' for more information.\n";
    return EXIT_FAILURE;
  }

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
  std::vector<std::pair<int, int>> su(1, std::make_pair(count_type(master, pebble::VertexType::MIN),
							count_type(master, pebble::VertexType::MAX)));
  std::vector<std::string> masters, reebs, morse_smales, duals;

  unsigned int index = 0;
  for (const MasterVertex vertex : boost::make_iterator_range(vertices(master)))
    master[vertex].id = index++;

  // walk the hierarchy of master graphs
  // and generate the Reeb and Morse-Smale graphs
  {
    const int min_count = count_type(master, pebble::VertexType::MIN);
    const int max_count = count_type(master, pebble::VertexType::MAX);
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

    const int min_count = count_type(master, pebble::VertexType::MIN);
    const int max_count = count_type(master, pebble::VertexType::MAX);

    su.emplace_back(min_count, max_count);
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
    difference.push_back(abs(s - stable_target) + abs(u - unstable_target));
  std::size_t min_index = std::distance(difference.cbegin(),
					min_element(difference.cbegin(),
					difference.cend()));

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
	    << masters[min_index] << ';'
	    << reebs[min_index] << ';'
	    << morse_smales[min_index] << ';'
	    << duals[min_index];

  return EXIT_SUCCESS;
}
