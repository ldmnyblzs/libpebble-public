#ifndef ENCODE_GRAPH_HPP
#define ENCODE_GRAPH_HPP 1

#include <bliss/graph.hh>
#include <sstream>
#include <set>

#include "types.hpp"

namespace pebble {
[[nodiscard]] std::string symbol(const unsigned int label) noexcept {
  if (label < 26) {
    return std::string(1, (char) ('a' + label));
  } else {
    return std::string(1, (char) ('A' + label / 26 - 1)).append(1, (char) ('a' + label % 26));
  }
}

template <typename Master>
std::string encode_graph(const Master &master) {
  using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
  using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;

  if (num_vertices(master) > 27 * 26)
    return std::string();

  bliss::Graph graph(num_vertices(master));

  for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
    graph.change_color(master[vertex].id, static_cast<unsigned int>(master[vertex].type));

  for (const MasterEdge &edge : boost::make_iterator_range(edges(master)))
    graph.add_edge(master[source(edge, master)].id,
		   master[target(edge, master)].id);

  bliss::Stats stats;
  const unsigned int *canonical = graph.canonical_form(stats, nullptr, nullptr);

  std::set<std::pair<unsigned int, unsigned int>> edge_set;

  for (const MasterEdge &edge : boost::make_iterator_range(edges(master)))
    edge_set.insert(std::minmax(canonical[master[source(edge, master)].id],
				canonical[master[target(edge, master)].id]));

  std::string result;
  for (const auto &[source, target] : edge_set)
    result.append(symbol(source)).append(symbol(target));

  return result;
}
}

#endif // ENCODE_GRAPH_HPP
