#ifndef CREATE_MASTER_GRAPH_HPP
#define CREATE_MASTER_GRAPH_HPP

#include <CGAL/boost/graph/properties.h>

#include "types.hpp"

namespace pebble {
template <typename Mesh,
	  typename PointMap,
	  typename EdgeMinimumMap,
	  typename FaceMinimumMap,
	  typename SaddleIterator,
	  typename AscendingIntersectionMap,
	  typename UnstableMap,
	  typename DescendingIntersectionMap,
	  typename StableMap,
	  typename Master,
	  typename VertexPointMap,
	  typename VertexTypeMap,
	  typename VertexSimplexMap,
	  typename EdgeTypeMap>
void create_master_graph(const Mesh &mesh,
			 const PointMap &point,
			 const EdgeMinimumMap &edge_minimum,
			 const FaceMinimumMap &face_minimum,
			 SaddleIterator saddle,
			 const SaddleIterator saddles_end,
			 const AscendingIntersectionMap &ascending_intersections,
			 const UnstableMap &unstables,
			 const DescendingIntersectionMap &descending_intersections,
			 const StableMap &stables,
			 Master &master,
			 VertexPointMap vertex_point,
			 VertexTypeMap vertex_type,
			 VertexSimplexMap vertex_simplex,
			 EdgeTypeMap edge_type) {
  using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
  using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
  using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
  using IntersectionSequence = typename boost::property_traits<AscendingIntersectionMap>::value_type;
  using Intersection = typename IntersectionSequence::value_type;

  std::map<MeshVertex, MasterVertex> vertex_to_vertex;
  std::map<MeshEdge, MasterVertex> edge_to_vertex;
  std::map<MeshFace, MasterVertex> face_to_vertex;

  for (; saddle != saddles_end; ++saddle) {
    const MasterVertex &saddle_vertex =
      get_or_add(edge_to_vertex, *saddle, master);
    put(vertex_point, saddle_vertex, get(edge_minimum, *saddle));
    put(vertex_type, saddle_vertex, VertexType::SADDLE);
    put(vertex_simplex, saddle_vertex, *saddle);

    const MeshHalfedge &h = halfedge(*saddle, mesh);
    for (const MeshHalfedge &halfedge : {h, opposite(h, mesh)}) {
      MasterVertex previous_vertex = saddle_vertex;
      for (const Intersection &up : get(ascending_intersections, halfedge)) {
	const MasterVertex new_intersection = add_vertex(master);
	put(vertex_point, new_intersection, up.at);
	put(vertex_type, new_intersection, VertexType::INTERSECTION);

	if (edge(new_intersection, previous_vertex, master).second)
	  std::cerr << "exists1\n";
	const MasterEdge ascending_edge =
	  add_edge(previous_vertex, new_intersection, master).first;
	put(edge_type, ascending_edge, EdgeType::ASCENDING);

	const MasterVertex &contour_saddle =
	  get_or_add(edge_to_vertex, edge(up.origin, mesh), master);
	if (edge(new_intersection, contour_saddle, master).second)
	  std::cerr << "exists2\n";
	const MasterEdge contour_edge =
	  add_edge(contour_saddle, new_intersection, master).first;
	put(edge_type, contour_edge, EdgeType::CONTOUR);

	previous_vertex = new_intersection;
      }
      const MeshVertex &unstable_mesh_vertex = get(unstables, halfedge);
      const MasterVertex &unstable_vertex =
	get_or_add(vertex_to_vertex, unstable_mesh_vertex, master);
      put(vertex_point, unstable_vertex, get(point, unstable_mesh_vertex));
      put(vertex_type, unstable_vertex, VertexType::MAX);
      put(vertex_simplex, unstable_vertex, unstable_mesh_vertex);

      if (edge(unstable_vertex, previous_vertex, master).second)
	std::cerr << "exists3\n";
      const MasterEdge ascending_edge =
	add_edge(previous_vertex, unstable_vertex, master).first;
      put(edge_type, ascending_edge, EdgeType::ASCENDING);

      const MeshFace &stable_mesh_face = get(stables, halfedge);
      previous_vertex = get_or_add(face_to_vertex, stable_mesh_face, master);
      put(vertex_point, previous_vertex, get(face_minimum, stable_mesh_face));
      put(vertex_type, previous_vertex, VertexType::MIN);
      put(vertex_simplex, previous_vertex, stable_mesh_face);
      for (const auto &down : get(descending_intersections, halfedge)) {
	const MasterVertex new_intersection = add_vertex(master);
	put(vertex_point, new_intersection, down.at);
	put(vertex_type, new_intersection, VertexType::INTERSECTION);

	if (edge(new_intersection, previous_vertex, master).second)
	  std::cerr << "exists4\n";
	const MasterEdge &ascending_edge =
	  add_edge(previous_vertex, new_intersection, master).first;
	put(edge_type, ascending_edge, EdgeType::ASCENDING);

	const MasterVertex &contour_saddle =
	  get_or_add(edge_to_vertex, edge(down.origin, mesh), master);
	if (edge(new_intersection, contour_saddle, master).second)
	  std::cerr << "exists5\n";
	const MasterEdge contour_edge =
	  add_edge(contour_saddle, new_intersection, master).first;
	put(edge_type, contour_edge, EdgeType::CONTOUR);

	previous_vertex = new_intersection;
      }
      auto exists = edge(saddle_vertex, previous_vertex, master);
      if (exists.second)
	std::cerr << "exists "
		  << ((int) get(vertex_type, source(exists.first, master))) << ' '
		  << ((int) get(vertex_type, target(exists.first, master))) << '\n';;
      const MasterEdge ascending_edge2 =
	add_edge(previous_vertex, saddle_vertex, master).first;
      put(edge_type, ascending_edge2, EdgeType::ASCENDING);
    }
  }
}
}

#endif // CREATE_MASTER_GRAPH_HPP
