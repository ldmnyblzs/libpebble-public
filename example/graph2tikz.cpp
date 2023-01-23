/*!
  \file graph2tikz.cpp
  \brief Generate a figure from an encoded string of a graph.
*/

#include <cstdlib>
#include <pebble/decode_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <sstream>
#include <fstream>

using Graph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			  //			  boost::property<boost::vertex_index_t, int>,
			  boost::no_property,
			  boost::property<boost::edge_index_t, int>>;
using GraphVertex = boost::graph_traits<Graph>::vertex_descriptor;
using GraphEdge = boost::graph_traits<Graph>::edge_descriptor;
using EdgeIterator = boost::graph_traits<Graph>::out_edge_iterator;
using Position = typename boost::convex_topology<2>::point_type;

/*!
  \brief Print the help message on the given output stream.
  \param stream the help text will be outputted to this stream
 */
void print_help(std::ostream &stream) {
  stream << "Usage: graph2tikz reeb|morse CODE\n"
	 << "Calculate the layout of the reeb or morse graph encoded as CODE.\n"
	 << "The resulting LaTeX file is printed to the standard output.\n";
}

int main(int argc, char *argv[]) {
  if (argc == 2 && std::strncmp(argv[1], "-h", 3) == 0) {
    print_help(std::cout);
    return EXIT_SUCCESS;
  }
  if (argc < 3) {
    print_help(std::cerr);
    return EXIT_FAILURE;
  }

  Graph graph;
  pebble::decode_graph(argv[2], graph);
  int index = 0;
  const auto edge_index = get(boost::edge_index, graph);
  for (const GraphEdge &edge : boost::make_iterator_range(edges(graph)))
    put(edge_index, edge, index++);
  const std::size_t n = num_vertices(graph);
  std::vector<Position> position_store(n);
  const auto position = boost::make_iterator_property_map(position_store.begin(),
							  get(boost::vertex_index, graph));

  if (strncmp(argv[1], "morse", 6) == 0) {
    boost::gursoy_atun_layout(graph,
			      boost::square_topology<>(1.3),
			      position);
    boost::fruchterman_reingold_force_directed_layout(graph,
    						      position,
    						      boost::square_topology<>(1.3));
  } else {
    const double radius = 1.0;
    boost::circle_graph_layout(graph, position, radius);
    boost::kamada_kawai_spring_layout(graph,
				      position,
				      boost::make_static_property_map<GraphEdge, double>(1.0),
				      boost::circle_topology<>(radius),
				      boost::edge_length(0.2));
  }

  std::cout << "\\documentclass[tikz]{standalone}\n"
	    << "\\begin{document}\n"
	    << "\\begin{tikzpicture}\n";

  for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
    const auto point = get(position, vertex);
    std::cout << "\\fill (" << point[0] << ',' << point[1] << ") circle (0.05);\n";
  }
  std::set<std::pair<GraphVertex, GraphVertex>> edge_visited;
  for (const GraphEdge &edge : boost::make_iterator_range(edges(graph))) {
    if (get(edge_index, edge) < index) {
      const GraphVertex s = source(edge, graph);
      const GraphVertex t = target(edge, graph);
      const auto point1 = get(position, s);
      const auto point2 = get(position, t);
      std::pair<GraphVertex, GraphVertex> e = std::minmax(s,t);
      const auto it = edge_visited.find(e);
      if (it == edge_visited.end() || *it != e) {
	edge_visited.emplace_hint(it, e);
	std::cout << "\\draw[gray] (" << point1[0] << ',' << point1[1]
		  << ") to (" << point2[0] << ',' << point2[1] << ");\n";
      } else {
	std::cout << "\\draw[gray, bend left] (" << point1[0] << ',' << point1[1]
		  << ") to (" << point2[0] << ',' << point2[1] << ");\n";
      }
    }
  }
  std::cout << "\\end{tikzpicture}\n"
	    << "\\end{document}\n";
  
  return EXIT_SUCCESS;
}
