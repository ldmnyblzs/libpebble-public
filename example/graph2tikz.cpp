/*!
  \file graph2tikz.cpp
  \brief Generate a figure from an encoded string of a graph.
*/

#include <boost/pending/property.hpp>
#include <cstdlib>
#include <pebble/decode_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/topology.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
// #include <boost/graph/random_layout.hpp>
// #include <boost/graph/fruchterman_reingold.hpp>

//#include "types.hpp"

#include <sstream>
#include <fstream>
#include <boost/graph/graphviz.hpp>

using Graph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			  //			  boost::property<boost::vertex_index_t, int>,
			  boost::no_property,
			  boost::property<boost::edge_index_t, int>>;
using GraphVertex = boost::graph_traits<Graph>::vertex_descriptor;
using GraphEdge = boost::graph_traits<Graph>::edge_descriptor;

using EmbeddingStorage = std::vector<std::vector<GraphEdge>>;
using IdMap = boost::property_map<Graph, boost::vertex_index_t>::type;
using EdgeIndexMap = boost::property_map<Graph, boost::edge_index_t>::type;
using Embedding =
    boost::iterator_property_map<EmbeddingStorage::iterator, IdMap>;

struct Point {
  std::size_t x;
  std::size_t y;
};
using DrawingStorage = std::vector<Point>;
using Drawing = boost::iterator_property_map<DrawingStorage::iterator, IdMap>;
using Position = typename boost::convex_topology<2>::point_type;

int main(int argc, char *argv[]) {
  if (argc < 3)
    return EXIT_FAILURE;

  Graph graph;
  pebble::decode_graph(argv[2], graph);
  int index = 0;
  EdgeIndexMap edge_index = get(boost::edge_index, graph);
  for (const GraphEdge &edge : boost::make_iterator_range(edges(graph)))
    put(edge_index, edge, index++);
  const std::size_t n = num_vertices(graph);
  std::vector<Position> position_store(n);
  auto position = boost::make_iterator_property_map(position_store.begin(),
						    get(boost::vertex_index, graph));

  if (strncmp(argv[1], "morse", 6) == 0) {
    Graph copy = graph;
    boost::edge_index_update_visitor<EdgeIndexMap> visitor(get(boost::edge_index, copy), index);
    boost::make_connected(copy, get(boost::vertex_index, copy), visitor);
    EmbeddingStorage embedding_storage(n);
    Embedding embedding(embedding_storage.begin(), get(boost::vertex_index, copy));
    if (!boost::boyer_myrvold_planarity_test(boost::boyer_myrvold_params::graph = copy,
					     boost::boyer_myrvold_params::embedding = embedding))
      return EXIT_FAILURE;
    boost::make_biconnected_planar(copy, embedding, get(boost::edge_index, copy), visitor);
    if (!boost::boyer_myrvold_planarity_test(boost::boyer_myrvold_params::graph = copy,
					     boost::boyer_myrvold_params::embedding = embedding))
      return EXIT_FAILURE;
    boost::make_maximal_planar(copy,
			       embedding,
			       get(boost::vertex_index, copy),
			       edge_index,
			       visitor);
    if (!boost::boyer_myrvold_planarity_test(boost::boyer_myrvold_params::graph = copy,
					     boost::boyer_myrvold_params::embedding = embedding))
      return EXIT_FAILURE;
    std::vector<GraphVertex> ordering;
    boost::planar_canonical_ordering(copy, embedding, std::back_inserter(ordering));
    DrawingStorage drawing_storage(n);
    Drawing drawing(drawing_storage.begin(), get(boost::vertex_index, copy));
    boost::chrobak_payne_straight_line_drawing(copy,
					       embedding,
					       ordering.begin(),
					       ordering.end(),
					       drawing);

    const double xmax = 2 * n - 4;
    const double ymax = n - 2;
    for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
      const Point &current = drawing[vertex];
      position[vertex][0] = current.x / xmax;
      position[vertex][1] = current.y / ymax;
    }
    
    // boost::circle_topology<> topology(1.0);
    // std::vector<double> weight_store(num_edges(graph));
    // auto weight = boost::make_iterator_property_map(weight_store.begin(),
    // 						    edge_index);
    // for (const GraphEdge &edge : boost::make_iterator_range(edges(graph)))
    //   put(weight, edge, get(edge_index, edge) < index ? 1.0 : sqrt(2.0));
    // boost::gursoy_atun_layout(graph,
    // 			      topology,
    // 			      position,
    // 			      1000 * num_vertices(graph),
    // 			      sqrt((double)num_vertices(graph)),
    // 			      1.0,
    // 			      0.8,
    // 			      0.2,
    // 			      get(boost::vertex_index, graph),
    // 			      weight);

    // for (int i = 0; i < num_vertices(graph); ++i) {
    //   for (int j = i + 1; j < num_vertices(graph); ++j) {
    // 	if (abs(position_store[i][0] - position_store[j][0]) < 0.01 &&
    // 	    abs(position_store[i][1] - position_store[j][1]) < 0.01) {
    // 	  position_store[i][0] += 0.1;
    // 	  position_store[i][1] += 0.1;
    // 	  break;
    // 	}
    //   }
    // }
    
    std::vector<double> weight_storage(num_edges(copy));
    auto weight = boost::make_iterator_property_map(weight_storage.begin(),
						    get(boost::edge_index, copy));
    GraphEdge last_edge;
    for (const GraphEdge &edge : boost::make_iterator_range(edges(copy))) {
      if (get(get(boost::edge_index, copy), edge) < index) {
	put(weight, edge, 1.0);
      } else {
	put(weight, edge, sqrt(2.0));
	last_edge = edge;
      }
    }
    remove_edge(last_edge, copy);
    const double side = 2.0;
    boost::square_topology<> topology(side);
    boost::kamada_kawai_spring_layout(copy,
				      position,
				      weight,
				      topology,
				      boost::side_length(side));
  } else {
    const double radius = 1.0;
    boost::circle_graph_layout(graph, position, radius);
    
    boost::circle_topology<> topology(radius);
    // boost::kamada_kawai_spring_layout(graph,
    // 				      position,
    // 				      boost::make_static_property_map<GraphEdge, double>(1.0),
    // 				      topology,
    // 				      boost::side_length(2 * radius),
    // 				      boost::layout_tolerance<double>(),
    // 				      10);
    boost::kamada_kawai_spring_layout(graph,
				      position,
				      boost::make_static_property_map<GraphEdge, double>(1.0),
				      topology,
				      boost::edge_length(0.2));
  }

  std::stringstream dotname;
  dotname << argv[2] << ".dot";
  std::ofstream dot(dotname.str());
  boost::write_graphviz(dot, graph);

  std::cout << "\\documentclass[tikz]{standalone}\n"
	    << "\\begin{document}\n"
	    << "\\begin{tikzpicture}\n";

  for (const GraphEdge &edge : boost::make_iterator_range(edges(graph))) {
    if (get(edge_index, edge) < index) {
      const auto point1 = get(position, source(edge, graph));
      const auto point2 = get(position, target(edge, graph));
      std::cout << "\\draw (" << point1[0] << ',' << point1[1]
		<< ") -- (" << point2[0] << ',' << point2[1] << ");\n";
    }
  }
  for (const GraphVertex &vertex : boost::make_iterator_range(vertices(graph))) {
    const auto point = get(position, vertex);
    std::cout << "\\fill (" << point[0] << ',' << point[1] << ") circle (0.1);\n";
  }
  std::cout << "\\end{tikzpicture}\n"
	    << "\\end{document}\n";
  
  return EXIT_SUCCESS;
}
