#ifndef DISTANCE_FROM_CENTROID_HPP
#define DISTANCE_FROM_CENTROID_HPP

#include <boost/graph/graph_traits.hpp>

namespace pebble {
template <typename Master,
	  typename Point,
	  typename PointMap,
	  typename DistanceMap>
void distance_from_centroid(const Master &master,
			    const Point &centroid,
			    const PointMap &point,
			    const DistanceMap &distance) {
  using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
  for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master))) {
    put(distance, vertex, sqrt(squared_distance(centroid, get(point, vertex))));
  }
}
}

#endif // DISTANCE_FROM_CENTROID_HPP
