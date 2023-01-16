/*!
  \file sort_saddles.hpp
  \brief Sort the saddle points by their persistence.
*/

#ifndef PEBBLE_SORT_SADDLES_HPP
#define PEBBLE_SORT_SADDLES_HPP 1

#include <boost/container/flat_map.hpp>

#include "types.hpp"
#include "utility.hpp"

namespace pebble {
  namespace internal {
    template <typename Master,
	      typename VertexTypeMap,
	      typename VertexDistanceMap,
	      typename EdgeTypeMap>
    typename boost::property_traits<VertexDistanceMap>::value_type
    descend(const Master &master,
	    const VertexTypeMap &vertex_type,
	    const VertexDistanceMap &distance,
	    const EdgeTypeMap &edge_type,
	    const typename boost::graph_traits<Master>::vertex_descriptor start) {
      using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
      using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
      using Scalar = typename boost::property_traits<VertexDistanceMap>::value_type;

      if (get(vertex_type, start) == VertexType::MIN)
	return get(distance, start);

      const MasterVertex saddle = source(in_edge_by_type(master, edge_type, start, EdgeType::CONTOUR), master);

      Scalar lowest = std::numeric_limits<Scalar>::max();
      for (const MasterEdge &edge : boost::make_iterator_range(in_edges(saddle, master))) {
	lowest = std::min(lowest,
			  descend(master,
				  vertex_type,
				  distance,
				  edge_type,
				  source(edge, master)));
      }
      return lowest;
    }

    template <typename Master, typename VertexTypeMap, typename VertexDistanceMap,
	      typename EdgeTypeMap>
    typename boost::property_traits<VertexDistanceMap>::value_type
    ascend(const Master &master,
	   const VertexTypeMap &vertex_type,
	   const VertexDistanceMap &distance,
	   const EdgeTypeMap &edge_type,
	   const typename boost::graph_traits<Master>::vertex_descriptor start) {
      using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
      using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
      using Scalar = typename boost::property_traits<VertexDistanceMap>::value_type;

      if (get(vertex_type, start) == VertexType::MAX)
	return get(distance, start);

      const MasterVertex saddle = source(in_edge_by_type(master, edge_type, start, EdgeType::CONTOUR), master);

      Scalar highest = std::numeric_limits<Scalar>::lowest();
      for (const MasterEdge &edge : boost::make_iterator_range(out_edges(saddle, master))) {
	if (get(edge_type, edge) == EdgeType::ASCENDING) {
	  highest = std::max(highest,
			     ascend(master,
				    vertex_type,
				    distance,
				    edge_type,
				    target(edge, master)));
	}
      }
      return highest;
    }
  };

  /*!
    \brief Sort the saddle points by their persistence.
    \tparam Master a model of BidirectionalGraph
    \tparam VertexTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and pebble::VertexType
    as value
    \tparam DistanceMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and a scalar with
    value
    \tparam EdgeTypeMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::edge_descriptor as key and pebble::EdgeType as
    value
    \tparam SaddleIterator a model of OutputIterator with
    boost::graph_traits<Master>::vertex_descriptor as value
    \param master the input graph
    \param vertex_type a property map containing the type of every vertex
    \param distance a property map containing the distance of the centroid and
    the equilibrium or intersection point vertices correspond to
    \param edge_type a property map containing the type of every edge
    \param saddle_out an iterator storing the saddles in ascending order of their persistence
  */
  template <typename Master,
	    typename VertexTypeMap,
	    typename DistanceMap,
	    typename EdgeTypeMap,
	    typename SaddleIterator>
  void sort_saddles(const Master &master,
		    const VertexTypeMap &vertex_type,
		    const DistanceMap &distance,
		    const EdgeTypeMap &edge_type,
		    SaddleIterator saddle_out) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
    using Scalar = typename boost::property_traits<DistanceMap>::value_type;
    using Map = boost::container::flat_map<MasterVertex, Scalar>;
    using Iterator = typename Map::iterator;
    using boost::container::small_vector;

    std::map<Scalar, MasterVertex> saddle_order;

    for (const MasterVertex &saddle : boost::make_iterator_range(vertices(master))) {
      if (get(vertex_type, saddle) == VertexType::SADDLE) {
	Map min_distance;
	for (const MasterEdge &edge : boost::make_iterator_range(in_edges(saddle, master))) {
	  const MasterVertex &adjacent = source(edge, master);
	  if (get(vertex_type, adjacent) == VertexType::MIN) { // stable
	    min_distance.emplace(adjacent, get(distance, adjacent));
	  } else { // intersection
	    const MasterVertex neighbor_saddle =
	      source(in_edge_by_type(master, edge_type, adjacent, EdgeType::CONTOUR), master);
	    Iterator iterator = min_distance.lower_bound(neighbor_saddle);
	    Scalar value = internal::descend(master, vertex_type, distance, edge_type, adjacent);
	    if (iterator == min_distance.end() || iterator->first != neighbor_saddle)
	      min_distance.emplace_hint(iterator, neighbor_saddle, value);
	    else
	      iterator->second = std::min(iterator->second, value);
	  }
	}

	Map max_distance;
	for (const MasterEdge &edge :
	       boost::make_iterator_range(out_edges(saddle, master))) {
	  if (get(edge_type, edge) == EdgeType::ASCENDING) {
	    const MasterVertex &adjacent = target(edge, master);
	    if (get(vertex_type, adjacent) == VertexType::MAX) { // unstable
	      max_distance.emplace(adjacent, get(distance, adjacent));
	    } else { // intersection
	      const MasterVertex neighbor_saddle =
		source(in_edge_by_type(master, edge_type, adjacent, EdgeType::CONTOUR), master);
	      Iterator iterator = max_distance.lower_bound(neighbor_saddle);
	      Scalar value = internal::ascend(master, vertex_type, distance, edge_type, adjacent);
	      if (iterator == max_distance.end() || iterator->first != neighbor_saddle)
		max_distance.emplace_hint(iterator, neighbor_saddle, value);
	      else
		iterator->second = std::max(iterator->second, value);
	    }
	  }
	}

	small_vector<Scalar, 3> distances;
	const Scalar &saddle_distance = get(distance, saddle);
	for (const auto [_, distance] : min_distance) {
	  distances.push_back(saddle_distance - distance);
	}
	for (const auto [_, distance] : max_distance) {
	  distances.push_back(distance - saddle_distance);
	}
	const Scalar min = *std::min_element(distances.cbegin(), distances.cend());
	saddle_order.emplace(min, saddle);
      }
    }

    for (const auto [distance, saddle] : saddle_order)
      *(saddle_out++) = saddle;
  }
};

#endif // PEBBLE_SORT_SADDLES_HPP
