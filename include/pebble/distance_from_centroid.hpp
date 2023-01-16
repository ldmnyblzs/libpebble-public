/*!
  \file distance_from_centroid.hpp
  \brief Calculate the distance of the centroid and the equilibrium or intersection point vertices correspond to in the master graph.
*/

#ifndef PEBBLE_DISTANCE_FROM_CENTROID_HPP
#define PEBBLE_DISTANCE_FROM_CENTROID_HPP 1

#include <boost/graph/graph_traits.hpp>

namespace pebble {
  /*!
    \brief Calculate the distance of the centroid and the equilibrium or
    intersection point vertices correspond to in the master graph.
    \tparam Master a model of VertexListGraph
    \tparam Point a model of Point_3
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and CGAL::Point_3 as
    value
    \tparam DistanceMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and a scalar with value
    \param master the input graph
    \param centroid the centroid of the mesh
    \param point a property map containing the position of the equilibrium or
    intersection point the vertex corresponds to
    \param distance a property map storing the distance of the centroid and the
    equilibrium or intersection point vertices correspond to
  */
  template <typename Master,
	    typename Point,
	    typename PointMap,
	    typename DistanceMap>
  void distance_from_centroid(const Master &master,
			      const Point &centroid,
			      const PointMap &point,
			      const DistanceMap &distance) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
      put(distance, vertex, sqrt(squared_distance(centroid, get(point, vertex))));
  }
};

#endif // PEBBLE_DISTANCE_FROM_CENTROID_HPP
