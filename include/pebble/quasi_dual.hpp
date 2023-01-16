/*!
  \file quasi_dual.hpp
  \brief Create the quasi dual corresponding to a Morse-Smale graph.
 */

#ifndef PEBBLE_QUASI_DUAL_HPP
#define PEBBLE_QUASI_DUAL_HPP 1

#include <unordered_map>
#include <boost/graph/graph_traits.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "types.hpp"

namespace pebble {
  /*!
    \brief Create the quasi dual corresponding to a Morse-Smale graph.
    \tparam Graph a model of MutableBidirectionalGraph
    \tparam VertexTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Graph>::vertex_descriptor as key and pebble::VertexType
    as value
    \param graph the input Morse-Smale graph
    \param vertex_type a property map containing the type of every vertex
    \returns the quasi dual corresponding to the Morse-Smale graph
  */
  template <typename Graph,
	    typename VertexTypeMap>
  Graph quasi_dual(const Graph &graph,
		   const VertexTypeMap &vertex_type) {
    using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
    using EdgeIterator = typename boost::graph_traits<Graph>::out_edge_iterator;

    Graph result;
    std::unordered_map<GraphVertex, GraphVertex> old_to_new;
    for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
      if (get(vertex_type, vertex) == VertexType::MAX) {
	const GraphVertex max = add_vertex(result);
	put(vertex_type, max, VertexType::MAX);
	old_to_new.emplace(vertex, max);
      }
    }
  
    for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
      if (get(vertex_type, vertex) == VertexType::MIN) {
	const GraphVertex new_min = add_vertex(result);
	put(vertex_type, new_min, VertexType::MIN);
	EdgeIterator edges_begin, edges_end;
	boost::tie(edges_begin, edges_end) = out_edges(vertex, graph);
	if (out_degree(vertex, graph) > 1) {
	  for (EdgeIterator edge1 = edges_begin;
	       boost::next(edge1) != edges_end;
	       edge1++) {
	    for (EdgeIterator edge2 = boost::next(edge1);
		 edge2 != edges_end;
		 edge2++) {
	      const GraphVertex saddle1 = target(*edge1, graph);
	      const GraphVertex saddle2 = target(*edge2, graph);
	      for (const GraphEdge &edge3 : boost::make_iterator_range(out_edges(saddle1, graph))) {
		for (const GraphEdge &edge4 : boost::make_iterator_range(out_edges(saddle2, graph))) {
		  const GraphVertex max = target(edge3, graph);
		  if (max == target(edge4, graph))
		    add_edge(new_min, old_to_new.find(max)->second, result);
		}
	      }
	    }
	  }
	} else if (out_degree(vertex, graph) == 1) {
	  const GraphVertex saddle = target(*edges_begin, graph);
	  const GraphEdge edge = *(out_edges(saddle, graph).first);
	  const GraphVertex max = target(edge, graph);
	  add_edge(new_min, old_to_new.find(max)->second, result);
	}
      }
    }
    return result;
  }
};

#endif // PEBBLE_QUASI_DUAL_HPP
