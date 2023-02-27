/*!
  \file pebble-cli.cpp
  \brief Classify a pebble stored in an STL file.
*/

#include <cstdlib>
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
#include <pebble/cancel_saddle3.hpp>
#include <pebble/create_reeb.hpp>
#include <pebble/encode_graph.hpp>
#include <pebble/create_morse_smale.hpp>
#include <pebble/get_saddles.hpp>
#include <pebble/create_master_graph_2.hpp>
#include <pebble/distance_from_centroid.hpp>
#include <pebble/axes.hpp>
#include <pebble/quasi_dual.hpp>
#include <pebble/measures.hpp>
//#include <pebble/classify_mesh_vertices.hpp>
//#include <pebble/classify_master_vertices.hpp>
#include <pebble/fit_paraboloid.hpp>

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
  bool convex;
  const bool read = pebble::convex_hull(argv[1], mesh, convex);
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
  const auto mesh_vertex_type =
    mesh.add_property_map<VertexHandle, pebble::VertexType>("v:type").first;
  const auto mesh_vertex_paramin =
    mesh.add_property_map<VertexHandle, Point>("v:paramin").first;

  // face and edge minimum
  pebble::face_center(mesh, centroid, face_center_pm, face_barycentric_pm, mesh.points());
  pebble::edge_saddle(mesh, centroid, edge_minimum_pm, edge_minimum_barycentric, mesh.points());

  // classify edges and vertices based on the position of the face minimum
  pebble::edge_type(mesh, face_barycentric_pm, is_followed_pm);
  pebble::vertex_dividing(mesh, face_barycentric_pm, dividing);

  // pebble::classify_mesh_vertices(mesh,
  // 				 centroid,
  // 				 get(boost::vertex_index, mesh),
  // 				 mesh_vertex_type,
  // 				 mesh_vertex_paramin);

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

  pebble::distance_from_centroid(master,
                                 centroid,
                                 get(&VertexData::point, master),
                                 get(&VertexData::distance, master));
  
  pebble::fit_paraboloid(mesh,
			 mesh.points(),
			 centroid,
			 master,
			 get(&VertexData::point, master),
			 get(&VertexData::type, master),
			 get(&VertexData::simplex, master),
			 get(&VertexData::keep, master));

  int stable = count_type(master, pebble::VertexType::MIN);
  int unstable = count_type(master, pebble::VertexType::MAX);
  int saddle = stable + unstable - 2;
  
  std::array<int, 3> type_count{0};
  for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
    if (master[vertex].type != pebble::VertexType::INTERSECTION &&
	master[vertex].keep)
      type_count[static_cast<int>(master[vertex].type)]++;
  const Scalar weight = static_cast<Scalar>(type_count[0]) / stable +
    static_cast<Scalar>(type_count[1]) / unstable +
    static_cast<Scalar>(type_count[2]) / saddle;
  
  std::map<Scalar, std::array<int, 2>> weights;
  weights.emplace(weight, std::array<int, 2>{stable, unstable});

  bool run = true;
  while (run && saddle > 2) {
    std::map<Scalar, std::array<MasterVertex, 3>> persistence;
    for (const MasterVertex &saddle : boost::make_iterator_range(vertices(master))) {
      if (master[saddle].type == pebble::VertexType::SADDLE) {
	std::vector<MasterVertex> stables;
	std::vector<MasterVertex> unstables;
	for (const MasterEdge &ascending : boost::make_iterator_range(out_edges(saddle, master))) {
	  if (master[ascending].type == pebble::EdgeType::ASCENDING) {
	    const MasterVertex adjacent = target(ascending, master);
	    switch (master[adjacent].type) {
	    case pebble::VertexType::MIN:
	      stables.push_back(adjacent);
	      break;
	    case pebble::VertexType::MAX:
	      unstables.push_back(adjacent);
	      break;
	    default:
	      break;
	    }
	  }
	}
	MasterVertex cancellable;
	MasterVertex kept;
	MasterVertex cancelled;
	bool found = false;
	if ((stables.size() == 2
	     && unstables.size() == 1
	     && stables.at(0) == stables.at(1)) ||
	    (stables.size() == 1
	     && unstables.size() == 2
	     && unstables.at(0) == unstables.at(1))) {
	  cancellable = saddle;
	  if ((master[saddle].distance - master[stables.at(0)].distance) >
	      (master[unstables.at(0)].distance - master[saddle].distance)) {
	    kept = stables.at(0);
	    cancelled = unstables.at(0);
	  } else {
	    kept = unstables.at(0);
	    cancelled = stables.at(0);
	  }
	  found = true;
	} else if (stables.size() == 2 && stables.at(0) != stables.at(1)) {
	  cancellable = saddle;
	  if (master[stables.at(0)].distance <
	      master[stables.at(1)].distance) {
	    kept = stables.at(0);
	    cancelled = stables.at(1);
	  } else {
	    kept = stables.at(1);
	    cancelled = stables.at(0);
	  }
	  found = true;
	} else if (unstables.size() == 2 && unstables.at(0) != unstables.at(1)) {
	  cancellable = saddle;
	  if (master[unstables.at(0)].distance >
	      master[unstables.at(1)].distance) {
	    kept = unstables.at(0);
	    cancelled = unstables.at(1);
	  } else {
	    kept = unstables.at(1);
	    cancelled = unstables.at(0);
	  }
	  found = true;
	}
	if (found)
	  persistence.emplace(abs(master[cancellable].distance - master[cancelled].distance),
			      std::array<MasterVertex, 3>{cancellable, kept, cancelled});
      }
    }
    const auto first = persistence.cbegin();
    if (first == persistence.cend()) {
      run = false;
    } else {
      const auto [cancellable, kept, cancelled] = first->second;
      pebble::cancel_saddle(master,
			    get(&VertexData::type, master),
			    get(&EdgeData::type, master),
			    cancellable,
			    kept,
			    cancelled);

      stable = count_type(master, pebble::VertexType::MIN);
      unstable = count_type(master, pebble::VertexType::MAX);
      saddle--;
      std::array<int, 3> type_count{0};
      for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
	if (master[vertex].type != pebble::VertexType::INTERSECTION &&
	    master[vertex].keep)
	  type_count[static_cast<int>(master[vertex].type)]++;
      const Scalar weight = static_cast<Scalar>(type_count[0]) / stable +
	static_cast<Scalar>(type_count[1]) / unstable +
	static_cast<Scalar>(type_count[2]) / saddle;
      
      weights.emplace(weight, std::array<int, 2>{stable, unstable});
    }
  }

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
	    << weights.rbegin()->second[0] << ';'
	    << weights.rbegin()->second[1] << ';'
	    << ";;;";
	    // << masters[min_index] << ';'
	    // << reebs[min_index] << ';'
	    // << morse_smales[min_index] << ';'
	    // << duals[min_index];

  return EXIT_SUCCESS;
}

	// if ((stables.size() == 2
	//      && unstables.size() == 1
	//      && stables.at(0) == stables.at(1)) ||
	//     (stables.size() == 1
	//      && unstables.size() == 2
	//      && unstables.at(0) == unstables.at(1))) {
	//   if (!master[saddle].keep) {
	//     cancellable = saddle;
	//     if (master[stables.at(0)].keep &&
	// 	!master[unstables.at(0)].keep) {
	//       kept = stables.at(0);
	//       cancelled = unstables.at(0);
	//     } else if (!master[stables.at(0)].keep &&
	// 	       master[unstables.at(0)].keep) {
	//       kept = unstables.at(0);
	//       cancelled = stables.at(0);
	//     } else if ((master[saddle].distance - master[stables.at(0)].distance) >
	// 	       (master[unstables.at(0)].distance - master[saddle].distance)) {
	//       kept = stables.at(0);
	//       cancelled = unstables.at(0);
	//     } else {
	//       kept = unstables.at(0);
	//       cancelled = stables.at(0);
	//     }
	//     run = true;
	//     //break;
	//   } else {
	//     if (master[stables.at(0)].keep &&
	// 	!master[unstables.at(0)].keep) {
	//       cancellable = saddle;
	//       kept = stables.at(0);
	//       cancelled = unstables.at(0);
	//       run = true;
	//       //break;
	//     } else if (!master[stables.at(0)].keep &&
	// 	       master[unstables.at(0)].keep) {
	//       cancellable = saddle;
	//       kept = unstables.at(0);
	//       cancelled = stables.at(0);
	//       run = true;
	//       //break;
	//     } else if (!master[stables.at(0)].keep &&
	// 	       !master[unstables.at(0)].keep) {
	//       cancellable = saddle;
	//       if ((master[saddle].distance - master[stables.at(0)].distance) >
	// 	  (master[unstables.at(0)].distance - master[saddle].distance)) {
	// 	kept = stables.at(0);
	// 	cancelled = unstables.at(0);
	//       } else {
	// 	kept = unstables.at(0);
	// 	cancelled = stables.at(0);
	//       }
	//       run = true;
	//       //break;
	//     }
	//   }
	// } else if (stables.size() == 2 && stables.at(0) != stables.at(1)) {
	//   if (!master[saddle].keep) {
	//     cancellable = saddle;
	//     if (master[stables.at(0)].keep &&
	// 	!master[stables.at(1)].keep) {
	//       kept = stables.at(0);
	//       cancelled = stables.at(1);
	//     } else if (!master[stables.at(0)].keep &&
	// 	       master[stables.at(1)].keep) {
	//       kept = stables.at(1);
	//       cancelled = stables.at(0);
	//     } else if (master[stables.at(0)].distance <
	// 	       master[stables.at(1)].distance) {
	//       kept = stables.at(0);
	//       cancelled = stables.at(1);
	//     } else {
	//       kept = stables.at(1);
	//       cancelled = stables.at(0);
	//     }
	//     run = true;
	//     //break;
	//   } else {
	//     if (master[stables.at(0)].keep &&
	// 	!master[stables.at(1)].keep) {
	//       cancellable = saddle;
	//       kept = stables.at(0);
	//       cancelled = stables.at(1);
	//       run = true;
	//       //break;
	//     } else if (!master[stables.at(0)].keep &&
	// 	       master[stables.at(1)].keep) {
	//       cancellable = saddle;
	//       kept = stables.at(1);
	//       cancelled = stables.at(0);
	//       run = true;
	//       //break;
	//     } else if (!master[stables.at(0)].keep &&
	// 	       !master[stables.at(1)].keep) {
	//       cancellable = saddle;
	//       if (master[stables.at(0)].distance <
	// 	  master[stables.at(1)].distance) {
	// 	kept = stables.at(0);
	// 	cancelled = stables.at(1);
	//       } else {
	// 	kept = stables.at(1);
	// 	cancelled = stables.at(0);
	//       }
	//       run = true;
	//       //break;
	//     }	    
	//   }
	// } else if (unstables.size() == 2 && unstables.at(0) != unstables.at(1)) {
	//   if (!master[saddle].keep) {
	//     cancellable = saddle;
	//     if (master[unstables.at(0)].keep &&
	// 	!master[unstables.at(1)].keep) {
	//       kept = unstables.at(0);
	//       cancelled = unstables.at(1);
	//     } else if (!master[unstables.at(0)].keep &&
	// 	       master[unstables.at(1)].keep) {
	//       kept = unstables.at(1);
	//       cancelled = unstables.at(0);
	//     } else if (master[unstables.at(0)].distance >
	// 	       master[unstables.at(1)].distance) {
	//       kept = unstables.at(0);
	//       cancelled = unstables.at(1);
	//     } else {
	//       kept = unstables.at(1);
	//       cancelled = unstables.at(0);
	//     }
	//     run = true;
	//     //break;
	//   } else {
	//     if (master[unstables.at(0)].keep &&
	// 	!master[unstables.at(1)].keep) {
	//       cancellable = saddle;
	//       kept = unstables.at(0);
	//       cancelled = unstables.at(1);
	//       run = true;
	//       //break;
	//     } else if (!master[unstables.at(0)].keep &&
	// 	       master[unstables.at(1)].keep) {
	//       cancellable = saddle;
	//       kept = unstables.at(1);
	//       cancelled = unstables.at(0);
	//       run = true;
	//       //break;
	//     } else if (!master[unstables.at(0)].keep &&
	// 	       !master[unstables.at(1)].keep) {
	//       cancellable = saddle;
	//       if (master[unstables.at(0)].distance >
	// 	  master[unstables.at(1)].distance) {
	// 	kept = unstables.at(0);
	// 	cancelled = unstables.at(1);
	//       } else {
	// 	kept = unstables.at(1);
	// 	cancelled = unstables.at(0);
	//       }
	//       run = true;
	//       //break;
	//     }
	//   }
	// }
