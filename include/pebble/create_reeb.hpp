/*!
  \file create_reeb.hpp
  \brief Create the Reeb-graph corresponding to a master graph.
*/

#ifndef PEBBLE_CREATE_REEB_HPP
#define PEBBLE_CREATE_REEB_HPP 1

#include "types.hpp"

namespace pebble {
  /*!
    \brief Create the Reeb-graph corresponding to a master graph.
    \tparam Master a model of MutableBidirectionalGraph
    \tparam EdgeTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::edge_descriptor as key and pebble::EdgeType as
    value
    \param master the input master graph
    \param edge_type a property map containing the type of every edge
    \return the Reeb-graph corresponding to the master graph
  */
  template <typename Master,
	    typename EdgeTypeMap>
  Master create_reeb(Master master,
		     const EdgeTypeMap &edge_type) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;

    std::vector<MasterVertex> to_delete;
    for (const MasterEdge &contour_edge : boost::make_iterator_range(edges(master))) {
      if (get(edge_type, contour_edge) == EdgeType::CONTOUR) {
	const MasterVertex deleted = target(contour_edge, master);
	const MasterVertex kept = source(contour_edge, master);
	for (const MasterEdge &in_edge : boost::make_iterator_range(in_edges(deleted, master))) {
	  const MasterVertex adjacent = source(in_edge, master);
	  if (adjacent != kept && !edge(adjacent, kept, master).second) {
	    const MasterEdge &edge = add_edge(adjacent, kept, master).first;
	    put(edge_type, edge, EdgeType::ASCENDING);
	  }
	}
	for (const MasterEdge &out_edge : boost::make_iterator_range(out_edges(deleted, master))) {
	  const MasterVertex adjacent = target(out_edge, master);
	  if (!edge(kept, adjacent, master).second) {
	    const MasterEdge &edge = add_edge(kept, adjacent, master).first;
	    put(edge_type, edge, EdgeType::ASCENDING);
	  }
	}
	to_delete.push_back(deleted);
      }
    }
    for (const MasterVertex &vertex : to_delete) {
      clear_vertex(vertex, master);
      remove_vertex(vertex, master);
    }
    return master;
  }
};

#endif // PEBBLE_CREATE_REEB_HPP
