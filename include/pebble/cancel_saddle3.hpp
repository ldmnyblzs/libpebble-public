#ifndef CANCEL_SADDLE_HPP
#define CANCEL_SADDLE_HPP 1

#include <ranges>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "types.hpp"
//#include "master_valid.hpp"

namespace pebble {
  template <typename Graph,
	    typename EdgeTypeMap>
  typename boost::graph_traits<Graph>::vertex_descriptor
  next_ascending(const Graph &graph,
		 const EdgeTypeMap &edge_type,
		 const typename boost::graph_traits<Graph>::vertex_descriptor &previous,
		 const typename boost::graph_traits<Graph>::vertex_descriptor &intersection) {
    using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

    for (const GraphEdge &ascending : make_iterator_range(out_edges(intersection, graph))) {
      if (get(edge_type, ascending) == EdgeType::ASCENDING) {
	const GraphVertex adjacent = opposite(ascending, intersection, graph);
	if (adjacent != previous) {
	  return adjacent;
	}
      }
    }
    throw "no ascending";
  }
  template <typename Graph,
	    typename VertexTypeMap,
	    typename EdgeTypeMap,
	    typename VertexIterator>
  void clear_intersections(Graph &graph,
			   const VertexTypeMap &vertex_type,
			   const EdgeTypeMap &edge_type,
			   const typename boost::graph_traits<Graph>::vertex_descriptor previous,
			   const typename boost::graph_traits<Graph>::vertex_descriptor current,
			   VertexIterator to_remove) {
    if (get(vertex_type, current) == pebble::VertexType::INTERSECTION) {
      *(to_remove++) = current;
      clear_intersections(graph,
			  vertex_type,
			  edge_type,
			  current,
			  next_ascending(graph, edge_type, previous, current),
			  to_remove);
    }
  }
  template <typename Graph,
	    typename VertexTypeMap,
	    typename EdgeTypeMap>
  void cancel_saddle(Graph &graph,
		     const VertexTypeMap &vertex_type,
		     const EdgeTypeMap &edge_type,
		     const typename boost::graph_traits<Graph>::vertex_descriptor saddle,
		     const typename boost::graph_traits<Graph>::vertex_descriptor kept,
		     const typename boost::graph_traits<Graph>::vertex_descriptor cancelled) {
    using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

    // BOOST_ASSERT(master_valid(graph, vertex_type, edge_type));
    // BOOST_ASSERT(cancelled != kept);
    // BOOST_ASSERT(edge(saddle, kept, graph).second);
    // BOOST_ASSERT(edge(saddle, cancelled, graph).second);
    // BOOST_ASSERT(get(edge_type, edge(saddle, kept, graph).first) == EdgeType::ASCENDING);
    // BOOST_ASSERT(get(edge_type, edge(saddle, cancelled, graph).first) == EdgeType::ASCENDING);

    // for (const GraphVertex &adjacent :
    // 	   boost::make_iterator_range(adjacent_vertices(cancelled, graph))) {
    //   if (adjacent != saddle) {
    // 	for (const GraphEdge &contour :
    // 	       boost::make_iterator_range(out_edges(adjacent, graph))) {
    // 	  if (get(edge_type, contour) == EdgeType::CONTOUR) {
    // 	    BOOST_ASSERT(boost::opposite(contour, adjacent, graph) == saddle);
    // 	  }
    // 	}
    //   }
    // }
    
    std::vector<GraphVertex> to_remove;
    for (const GraphEdge &edge :
	   boost::make_iterator_range(out_edges(saddle, graph))) {
      const GraphVertex adjacent = opposite(edge, saddle, graph);
      switch (get(edge_type, edge)) {
      case EdgeType::ASCENDING:
	if (adjacent != kept && adjacent != cancelled) {
	  clear_intersections(graph,
			      vertex_type,
			      edge_type,
			      saddle,
			      adjacent,
			      std::back_inserter(to_remove));
	}
	break;
      case EdgeType::CONTOUR: {
	for (const GraphEdge &ascending : make_iterator_range(out_edges(adjacent, graph))) {
	  if (get(edge_type, ascending) == EdgeType::ASCENDING) {
	    const GraphVertex two_ring = opposite(ascending, adjacent, graph);
	    if (two_ring != kept && two_ring != cancelled) {
	      const std::pair<GraphEdge, bool> newedge = add_edge(two_ring, kept, graph);
	      put(edge_type,
		  newedge.first,
		  EdgeType::ASCENDING);
	      to_remove.push_back(adjacent);
	    }
	  }
	}
      }
	break;
      }
    }

    to_remove.push_back(saddle);
    to_remove.push_back(cancelled);
    std::sort(to_remove.begin(),
	      to_remove.end(),
	      std::greater<GraphVertex>());
    for (const GraphVertex &vertex : to_remove) {
      clear_vertex(vertex, graph);
    }
    for (const GraphVertex &vertex : to_remove) {
      remove_vertex(vertex, graph);
    }
  }
};

#endif // CANCEL_SADDLE_HPP
