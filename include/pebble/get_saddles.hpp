/*!
  \file get_saddles.hpp
  \brief Get all saddles on the input mesh.
*/

#ifndef PEBBLE_GET_SADDLES
#define PEBBLE_GET_SADDLES

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

namespace pebble {
  /*!
    \brief Get all saddles on the input mesh.
    \tparam Mesh model of EdgeListGraph
    \tparam EdgeFollowedPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and bool as value
    \tparam EdgeMinimumBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and a Container of 2 scalar
    values as value
    \tparam SaddleIterator a model of OutputIterator with
    boost::graph_traits<Mesh>::edge_descriptor as value
    \param mesh the input mesh
    \param edge_followed a property map containing true as value for followed edges and
    false for crossed edges
    \param edge_minimum_barycentric a property map containing the position of the centroid projected to the line of each edge in barycentric coordinates
    \param saddles iterator taking the edge descriptors of edges with saddles
  */
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
};

#endif // PEBBLE_GET_SADDLES
