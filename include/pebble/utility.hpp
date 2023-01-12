#ifndef UTILITY_HPP
#define UTILITY_HPP 1

#include "types.hpp"
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/property_map/property_map.hpp>

namespace pebble {

template <typename Kernel>
[[nodiscard]] bool opposite_direction(const CGAL::Vector_3<Kernel> &v1,
				      const CGAL::Vector_3<Kernel> &v2) noexcept {
  return std::make_tuple(std::signbit(v1.x()),
			 std::signbit(v1.y()),
			 std::signbit(v1.z())) !=
    std::make_tuple(std::signbit(v2.x()),
		    std::signbit(v2.y()),
		    std::signbit(v2.z()));
}

template <typename Kernel>
[[nodiscard]] typename Kernel::FT angle_on_plane(const CGAL::Point_3<Kernel> &center,
						 const CGAL::Vector_3<Kernel> &normal,
						 const CGAL::Point_3<Kernel> &from,
						 const CGAL::Point_3<Kernel> &to) noexcept {
  //std::cout << center << '\t' << from << '\t' << to << std::endl;
  if (from == to || from == center || to == center)
    return 360;
  const typename Kernel::FT angle = CGAL::approximate_angle(from, center, to);
  if (opposite_direction(normal, CGAL::unit_normal(center, from, to)))
    return 360 - angle;
  return angle;
}

template <typename Mesh, typename Kernel>
[[nodiscard]] bool
is_ascending(const Mesh &mesh,
	     const CGAL::Point_3<Kernel> &centroid,
	     const typename boost::graph_traits<Mesh>::halfedge_descriptor &halfedge) noexcept {
  return CGAL::has_smaller_distance_to_point(centroid,
					     mesh.point(source(halfedge, mesh)),
					     mesh.point(target(halfedge, mesh)));
}

template <typename Kernel>
[[nodiscard]] typename Kernel::FT
arc_length(const CGAL::Point_3<Kernel> &center,
	   const CGAL::Vector_3<Kernel> &normal,
	   const CGAL::Point_3<Kernel> &from,
	   const CGAL::Point_3<Kernel> &to) noexcept {
  return angle_on_plane(center, normal, from, to) * sqrt(squared_distance(center, from));
}

template <typename Kernel>
[[nodiscard]] CGAL::Vector_3<Kernel>
rotate(const CGAL::Vector_3<Kernel> &vector,
       const CGAL::Vector_3<Kernel> &normal,
       const typename Kernel::FT &angle) noexcept {
  using Scalar = typename Kernel::FT;

  const Scalar radians = angle * boost::math::constants::pi<Scalar>() / 180.0;
  const Scalar sine = sin(radians);
  const Scalar cosine = cos(radians);
  return vector * cosine
    + cross_product(normal, vector) * sine
    + normal * scalar_product(normal, vector) * (1 - cosine);
}

template <typename Scalar>
Scalar to_radians(const Scalar &degrees) {
  return degrees / static_cast<Scalar>(180.0) * boost::math::constants::pi<Scalar>();
}

template <typename Vector>
Vector normalize(const Vector &vector) {
  return vector / sqrt(vector.squared_length());
}

template <typename Map, typename Master>
typename boost::graph_traits<Master>::vertex_descriptor
get_or_add(Map &map,
	   const typename Map::key_type &key,
	   Master &master) {
  const auto it = map.lower_bound(key);
  if (it != map.end() && it->first == key)
    return it->second;

  return map.emplace_hint(it,
			  key,
			  add_vertex(master))->second;
}

template <typename Graph,
	  typename EdgeTypeMap,
	  typename EdgeUnderlyingMap = boost::typed_identity_property_map<typename boost::graph_traits<Graph>::edge_descriptor>>
typename boost::graph_traits<Graph>::edge_descriptor
in_edge_by_type(const Graph &graph,
		const EdgeTypeMap &edge_type,
		const typename boost::graph_traits<Graph>::vertex_descriptor &vertex,
		const EdgeType type,
		const EdgeUnderlyingMap &underlying = boost::typed_identity_property_map<typename boost::graph_traits<Graph>::edge_descriptor>()) {
  using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
  for (const GraphEdge &edge : boost::make_iterator_range(in_edges(vertex, graph)))
    if (get(edge_type, get(underlying, edge)) == type)
      return edge;
  throw "error";
}

template <typename Graph,
	  typename EdgeTypeMap,
	  typename EdgeUnderlyingMap = boost::typed_identity_property_map<typename boost::graph_traits<Graph>::edge_descriptor>>
typename boost::graph_traits<Graph>::edge_descriptor
out_edge_by_type(const Graph &graph,
		 const EdgeTypeMap &edge_type,
		 const typename boost::graph_traits<Graph>::vertex_descriptor &vertex,
		 const EdgeType type,
		 const EdgeUnderlyingMap &underlying = boost::typed_identity_property_map<typename boost::graph_traits<Graph>::edge_descriptor>()) {
  using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
  for (const GraphEdge &edge : boost::make_iterator_range(out_edges(vertex, graph)))
    if (get(edge_type, get(underlying, edge)) == type)
      return edge;
  throw "error";
}
}

#endif // UTILITY_HPP
