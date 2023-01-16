/*!
  \file create_master_graph.hpp
  \brief Create the master graph based on the equilibrium and intersection
  points on the surface of the mesh.
*/

#ifndef PEBBLE_CREATE_MASTER_GRAPH_HPP
#define PEBBLE_CREATE_MASTER_GRAPH_HPP

#include <CGAL/boost/graph/properties.h>

#include "types.hpp"

namespace pebble {
  /*!
    \brief Create the master graph based on the equilibrium and intersection
    points on the surface of the mesh.
    \tparam Mesh a model of FaceGraph
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as
    value
    \tparam EdgeMinimumMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and CGAL::Point_3 as value
    \tparam FaceMinimumMap a model of ReadablePropertyMap
    with boost::graph_traits<Mesh>::face_descriptor as key and CGAL::Point_3 as
    value
    \tparam SaddleIterator a model of InputIterator with
    boost::graph_traits<Mesh>::edge_descriptor as value
    \tparam IntersectionMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and a vector of
    pebble::Intersection objects as value
    \tparam UnstableMap a model of ReadablePropertyMap
    with boost::graph_traits<Mesh>::halfedge_descriptor as key and
    boost::graph_traits<Mesh>::vertex_descriptor as value
    \tparam StableMap a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and
    boost::graph_traits<Mesh>::face_descriptor as value
    \tparam Master a model of MutableGraph
    \tparam VertexPointMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and CGAL::Point_3 as
    value
    \tparam VertexTypeMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and pebble::VertexType
    as value
    \tparam VertexSimplexMap a model of WritablePropertyMap
    boost::graph_traits<Master>::vertex_descriptor as key and
    std::variant<std::monostate, boost::graph_traits<Mesh>::vertex_descriptor,
    boost::graph_traits<Mesh>::edge_descriptor,
    boost::graph_traits<Mesh>::face_descriptor> as value
    \tparam EdgeTypeMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::edge_descriptor as key and pebble::EdgeType as
    value
    \param mesh the input mesh
    \param point a property map providing the position of each vertex of the
    mesh
    \param edge_minimum  a property map containing the position of the centroid
    projected to the line of each edge of the mesh in Descartes coordinates
    \param face_minimum a property map containing the position of the
    centroid projected to the plane of each face of the mesh in Descartes
    coordinates
    \param saddle the beginning of the range of edges with saddle
    \param saddles_end the end of the range of edges with saddle
    \param ascending_intersections a property map containing for every halfedge
    the intersection points on the saddle-unstable isolated ascending curve
    originating from the saddle on the halfedge
    \param unstables a property map containing for every halfedge the unstable
    point at the end of the isolated ascending curve originating from the saddle
    on the halfedge and ascending in the direction of the halfedge's target
    \param descending_intersections a property map containing for every halfedge
    the intersection points on the stable-saddle isolated ascending curve
    originating from the saddle on the halfedge
    \param stables  a property map containing for every halfedge the stable
    point at the end of the isolated ascending curve originating from the saddle
    on the halfedge and descending in the direction of the face adjacent to the
    halfedge
    \param master the output graph
    \param vertex_point a property map storing the position of the equilibrium
    or intersection point the vertex corresponds to
    \param vertex_type a property map storing the type of every vertex
    \param vertex_simplex a property map storing the vertex, edge or face of the
    mesh that contains the equilibrium or intersection point the vertex corresponds
    to
    \param edge_type a property map storing the type of every edge
  */
  template <typename Mesh,
	    typename PointMap,
	    typename EdgeMinimumMap,
	    typename FaceMinimumMap,
	    typename SaddleIterator,
	    typename IntersectionMap,
	    typename UnstableMap,
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
			   const IntersectionMap &ascending_intersections,
			   const UnstableMap &unstables,
			   const IntersectionMap &descending_intersections,
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
    using IntersectionSequence = typename boost::property_traits<IntersectionMap>::value_type;
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

	  const MasterEdge ascending_edge =
	    add_edge(previous_vertex, new_intersection, master).first;
	  put(edge_type, ascending_edge, EdgeType::ASCENDING);

	  const MasterVertex &contour_saddle =
	    get_or_add(edge_to_vertex, edge(up.origin, mesh), master);
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

	  const MasterEdge &ascending_edge =
	    add_edge(previous_vertex, new_intersection, master).first;
	  put(edge_type, ascending_edge, EdgeType::ASCENDING);

	  const MasterVertex &contour_saddle =
	    get_or_add(edge_to_vertex, edge(down.origin, mesh), master);
	  const MasterEdge contour_edge =
	    add_edge(contour_saddle, new_intersection, master).first;
	  put(edge_type, contour_edge, EdgeType::CONTOUR);

	  previous_vertex = new_intersection;
	}
	auto exists = edge(saddle_vertex, previous_vertex, master);
	const MasterEdge ascending_edge2 =
	  add_edge(previous_vertex, saddle_vertex, master).first;
	put(edge_type, ascending_edge2, EdgeType::ASCENDING);
      }
    }
  }
};

#endif // PEBBLE_CREATE_MASTER_GRAPH_HPP
