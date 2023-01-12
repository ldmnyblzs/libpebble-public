#ifndef PEBBLE_DECODE_GRAPH_HPP
#define PEBBLE_DECODE_GRAPH_HPP 1

#include <string>
#include <vector>
#include <array>
#include <optional>
#include <boost/graph/graph_traits.hpp>

namespace pebble {
template <typename Master>
void decode_graph(const std::string& encoded,
		  Master &master) {
  using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
  using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
  
  std::vector<std::array<unsigned int, 2>> edges;

  std::optional<unsigned int> previous;
  const std::size_t length = encoded.size();
  std::size_t i = 0;
  unsigned int last_vertex = 0;
  while(i < length) {
    unsigned int vertex = 0;
    if (encoded[i] <= 'Z') {
      vertex += (encoded[i] - 'A' + 1) * 26 + (encoded[i + 1] - 'a');
      i += 2;
    } else {
      vertex += encoded[i] - 'a';
      ++i;
    }
    last_vertex = std::max(last_vertex, vertex);
    if (previous.has_value()) {
      edges.push_back({previous.value(), vertex});
      previous.reset();
    } else {
      previous = vertex;
    }
  }

  std::vector<MasterVertex> vertices;
  for (unsigned int i = 0; i <= last_vertex; ++i) {
    const MasterVertex vertex = add_vertex(master);
    vertices.push_back(vertex);
    //put(index, vertex, i);
  }

  for (const auto &[source, target] : edges)
    add_edge(vertices.at(source), vertices.at(target), master);
}
};

#endif // PEBBLE_DECODE_GRAPH_HPP
