#ifndef GET_SADDLES
#define GET_SADDLES

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

namespace pebble {
template <typename Mesh,
	  typename EdgeFollowedPM,
	  typename EdgeMinimumBarycentricPM,
	  typename SaddleIterator>
void get_saddles(const Mesh &mesh,
		 const EdgeFollowedPM &edge_followed,
		 const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		 SaddleIterator saddles) {
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using Barycentric = typename boost::property_traits<EdgeMinimumBarycentricPM>::value_type;

  for (const MeshEdge &edge : edges(mesh)) {
    const Barycentric &minimum = get(edge_minimum_barycentric, edge);
    if (get(edge_followed, edge) &&
	!(minimum[0] < 0 ||
	 minimum[1] < 0)) {
      *(saddles++) = edge;
    }
  }
}
}

#endif // GET_SADDLES
