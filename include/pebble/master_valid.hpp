#ifndef MASTER_VALID_HPP
#define MASTER_VALID_HPP 1

#include <boost/graph/graph_traits.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <map>
#include <optional>

#include "types.hpp"

namespace pebble {
  template <typename Graph,
	    typename EdgeTypeMap>
  typename boost::graph_traits<Graph>::vertex_descriptor
  contour(const Graph &graph,
	  const EdgeTypeMap &edge_type,
	  const typename boost::graph_traits<Graph>::vertex_descriptor &intersection) {
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
    for (const GraphEdge &contour : make_iterator_range(out_edges(intersection, graph))) {
      if (get(edge_type, contour) == EdgeType::CONTOUR) {
	return opposite(contour, intersection, graph);
      }
    }
    throw "no contour";
  }
  template <typename Graph,
	    typename VertexTypeMap,
	    typename EdgeTypeMap>
  typename boost::graph_traits<Graph>::vertex_descriptor
  ascend(const Graph &graph,
	 const VertexTypeMap &vertex_type,
	 const EdgeTypeMap &edge_type,
	 const typename boost::graph_traits<Graph>::vertex_descriptor &previous,
	 const typename boost::graph_traits<Graph>::vertex_descriptor &current) {
    using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
    if (get(vertex_type, current) != VertexType::INTERSECTION) {
      return current;
    }
    for (const GraphEdge &ascending : make_iterator_range(out_edges(current, graph))) {
      if (get(edge_type, ascending) == EdgeType::ASCENDING) {
	const GraphVertex adjacent = opposite(ascending, current, graph);
	if (adjacent != previous) {
	  return ascend(graph,
			vertex_type,
			edge_type,
			current,
			adjacent);
	}
      }
    }
    std::cout << graph[current].id << ' ';
    throw "no ascending";
  }
  
  template <typename Graph,
	    typename VertexTypeMap,
	    typename EdgeTypeMap>
  bool
  master_valid(const Graph &graph,
	       const VertexTypeMap &vertex_type,
	       const EdgeTypeMap &edge_type) {    
    using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;
    unsigned int zero = 0;
    bool result = true;
    for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
      if (graph[vertex].id == 0)
	++zero;
      unsigned int a = 0;
      unsigned int c = 0;
      for (const GraphEdge &edge : make_iterator_range(out_edges(vertex, graph))) {
	switch (get(edge_type, edge)) {
	case EdgeType::ASCENDING: {
	  std::array<GraphVertex, 2> v = {source(edge, graph), opposite(edge, source(edge, graph), graph)};
	  try {
	    if (get(vertex_type, v[0]) == VertexType::INTERSECTION &&
		get(vertex_type, v[1]) == VertexType::INTERSECTION &&
		contour(graph, edge_type, v[0]) == contour(graph, edge_type, v[1])) {
	      std::cout << 'a' << ' '
			<< contour(graph, edge_type, v[0]) << ' '
			<< contour(graph, edge_type, v[1]) << std::endl;
	      result = false;
	    }
	  } catch (const char *error) {
	    std::cout << error << std::endl;
	  }
	  ++a;
	  break;
	}
	case EdgeType::CONTOUR:
	  ++c;
	  break;
	default:
	  std::cout << 'x';
	  result = false;
	}
      }
      switch (get(vertex_type, vertex)) {
      case VertexType::INTERSECTION:
	if (a != 2) {
	  std::cout << 'b' << a << std::endl;
	  result = false;
	}
	if (c != 1) {
	  std::cout << 'c' << c << std::endl;
	  result = false;
	}
	try {
	  std::optional<GraphVertex> other;
	  for (const GraphEdge &ascending :
		 boost::make_iterator_range(out_edges(vertex, graph))) {
	    if (get(edge_type, ascending) == EdgeType::ASCENDING) {
	      const GraphVertex adjacent = opposite(ascending, vertex, graph);
	      const GraphVertex elliptic = ascend(graph,
						  vertex_type,
						  edge_type,
						  vertex,
						  adjacent);
	      if (other.has_value()) {
		if (get(vertex_type, *other) == get(vertex_type, elliptic)) {
		  std::cout << 'g' << std::endl;
		  result = false;
		}
	      } else {
		other = elliptic;
	      }
	    }
	  }
	  if (a == 1)
	    std::cout << get(vertex_type, *other) << *other << std::endl;
	} catch (const char *error) {
	  
	}
	break;
      case VertexType::MIN:
      case VertexType::MAX: {
	if (a == 0) {
	  std::cout << 'd' << c << std::endl;
	  result = false;
	}
	if (c != 0) {
	  std::cout << 'e' << c << std::endl;
	  result = false;
	}
	std::optional<GraphVertex> crossed;
	for (const GraphVertex &adjacent :
	       boost::make_iterator_range(adjacent_vertices(vertex, graph))) {
	  const GraphVertex saddle = ascend(graph,
					    vertex_type,
					    edge_type,
					    vertex,
					    adjacent);
	  if (get(vertex_type, saddle) != VertexType::SADDLE) {
	    std::cout << 'h' << get(vertex_type, saddle) << std::endl;
	    result = false;
	  }
	  switch (get(vertex_type, adjacent)) {
	  case VertexType::INTERSECTION: {
	    const GraphVertex cont = contour(graph, edge_type, adjacent);
	    if (crossed.has_value()) {
	      if (crossed != cont) {
		std::cout << 'k'
			  << get(vertex_type, vertex) << graph[vertex].id
			  << '-' << graph[*crossed].id
			  << '-' << graph[cont].id << std::endl;
		result = false;
	      }
	    } else {
	      crossed = cont;
	    }
	    break;
	  }
	  case VertexType::SADDLE: {
	    if (crossed.has_value()) {
	      if (crossed != adjacent) {
		std::cout << 'l' << std::endl;
		result = false;
	      }
	    } else {
	      crossed = adjacent;
	    }
	    break;
	  }
	  }
	}
	break;
	}
      case VertexType::SADDLE:
	if (a != 4) {
	  std::cout << 'f' << a << ' ' << graph[vertex].id << std::endl;
	  result = false;
	}
	try {
	  unsigned short S = 0, U = 0;
	  for (const GraphEdge &ascending :
		 boost::make_iterator_range(out_edges(vertex, graph))) {
	    if (get(edge_type, ascending) == EdgeType::ASCENDING) {
	      const GraphVertex adjacent = opposite(ascending, vertex, graph);
	      const GraphVertex elliptic = ascend(graph,
						  vertex_type,
						  edge_type,
						  vertex,
						  adjacent);
	      switch (get(vertex_type, elliptic)) {
	      case VertexType::MIN:
		++S;
		break;
	      case VertexType::MAX:
		++U;
		break;
	      default:
		std::cout << 'i' << std::endl;
		result = false;
	      }
	    }
	  }
	  if (S != 2 || U != 2) {
	    std::cout << 'j' << S << ' ' << U << std::endl;
	    result = false;
	  }
	} catch (const char *error) {
	  std::cout << error << std::endl;
	}
	break;
      default:
	std::cout << 'X';
      }
    }
    if (zero > 1) {
      std::cout << 'z' << std::endl;
      result = false;
    }
    return result;
  }
};

#endif // MASTER_VALID_HPP
