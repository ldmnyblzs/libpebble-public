/*!
  \file create_morse_smale.hpp
  \brief Create the Morse-Smale graph corresponding to a master graph.
*/

#ifndef PEBBLE_CREATA_MORSE_SMALE_HPP
#define PEBBLE_CREATA_MORSE_SMALE_HPP 1

#include "types.hpp"

namespace pebble {
  namespace internal {
    //! Clear an intersection and connect the vertices of the two adjacent ascending edges opposite to the intersection.
    template <typename Master,
	      typename EdgeTypeMap>
    void clear_intersection(Master &master,
			    const typename boost::graph_traits<Master>::vertex_descriptor &intersection,
			    const EdgeTypeMap &edge_type) {
      using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
      using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;

      const MasterEdge &out_edge = *out_edges(intersection, master).first;
      const MasterVertex &t = target(out_edge, master);
      MasterEdge in_edge;
      for (const MasterEdge &edge : boost::make_iterator_range(in_edges(intersection, master)))
	if (get(edge_type, edge) == EdgeType::ASCENDING)
	  in_edge = edge;
      const MasterVertex &s = source(in_edge, master);
      const MasterEdge &edge = add_edge(s, t, master).first;
      put(edge_type, edge, EdgeType::ASCENDING);
      clear_vertex(intersection, master);
    }
  };

  /*!
    \brief Create the Morse-Smale graph corresponding to a master graph.
    \tparam Master a model of MutableBidirectionalGraph
    \tparam VertexTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and pebble::VertexType
    as value
    \tparam EdgeTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::edge_descriptor as key and pebble::EdgeType as
    value
    \param master the input master graph
    \param vertex_type a property map containing the type of every vertex
    \param edge_type a property map containing the type of every edge
    \return the Morse-Smale graph corresponding to the master graph
   */
  template <typename Master,
	    typename VertexTypeMap,
	    typename EdgeTypeMap>
  Master create_morse_smale(Master master,
			    const VertexTypeMap &vertex_type,
			    const EdgeTypeMap &edge_type) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;

    std::vector<MasterVertex> isolated;
    for (const MasterVertex &intersection : boost::make_iterator_range(vertices(master))) {
      if (get(vertex_type, intersection) == VertexType::INTERSECTION) {
	internal::clear_intersection(master, intersection, edge_type);
	isolated.push_back(intersection);
      }
    }
    for (const MasterVertex &vertex : isolated)
      remove_vertex(vertex, master);
    return master;
  }
};

#endif // PEBBLE_CREATA_MORSE_SMALE_HPP
