/*!
  \file cancel_saddle.hpp
  \brief Cancel a saddle vertex of the master graph.
*/

#ifndef PEBBLE_CANCEL_SADDLE_HPP
#define PEBBLE_CANCEL_SADDLE_HPP
#include "boost/graph/properties.hpp"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/reverse_graph.hpp>

#include "utility.hpp"

namespace pebble {
  namespace internal {
    template <typename Graph,
	      typename VertexTypeMap,
	      typename EdgeUnderlyingMap,
	      typename EdgeTypeMap,
	      typename IntersectionIterator>
    typename boost::graph_traits<Graph>::vertex_descriptor
    intersections_up(const Graph &graph,
		     const VertexTypeMap &vertex_type,
		     const EdgeUnderlyingMap &underlying,
		     const EdgeTypeMap &edge_type,
		     const typename boost::graph_traits<Graph>::vertex_descriptor &from,
		     IntersectionIterator intersection) {
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      if (get(vertex_type, from) == VertexType::INTERSECTION) {
	*(intersection++) = from;
	const GraphEdge &ascending = out_edge_by_type(graph,
						      edge_type,
						      from,
						      EdgeType::ASCENDING,
						      underlying);
	return intersections_up(graph,
				vertex_type,
				underlying,
				edge_type,
				target(ascending, graph),
				intersection);
      }
      return from;
    }

    template <typename Graph,
	      typename VertexTypeMap,
	      typename EdgeUnderlyingMap,
	      typename EdgeTypeMap,
	      typename IntersectionIterator>
    void
    intersections_down(const Graph &graph,
		       const VertexTypeMap &vertex_type,
		       const EdgeUnderlyingMap &underlying,
		       const EdgeTypeMap &edge_type,
		       const typename boost::graph_traits<Graph>::vertex_descriptor &from,
		       IntersectionIterator intersection) {
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      if (get(vertex_type, from) == VertexType::INTERSECTION) {
	*(intersection++) = from;
	const GraphEdge &ascending = in_edge_by_type(graph,
						     edge_type,
						     from,
						     EdgeType::ASCENDING,
						     underlying);
	intersections_down(graph,
			   vertex_type,
			   underlying,
			   edge_type,
			   source(ascending, graph),
			   intersection);
      }
    }

    template<typename Graph,
	     typename VertexTypeMap,
	     typename VertexDistanceMap,
	     typename EdgeUnderlyingMap,
	     typename EdgeTypeMap,
	     typename Compare,
	     typename IntersectionIterator>
    std::pair<typename boost::graph_traits<Graph>::vertex_descriptor,
	      typename boost::graph_traits<Graph>::vertex_descriptor>
    find_unstables(const Graph &graph,
		   const VertexTypeMap &vertex_type,
		   const VertexDistanceMap &distance,
		   const EdgeUnderlyingMap &underlying,
		   const EdgeTypeMap &edge_type,
		   const typename boost::graph_traits<Graph>::vertex_descriptor &saddle,
		   const VertexType unstable_type,
		   Compare compare,
		   IntersectionIterator intersection) {
      using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      std::optional<GraphVertex> lower;
      GraphVertex higher;
      for (const GraphEdge &edge : boost::make_iterator_range(out_edges(saddle, graph))) {
	if (get(edge_type, get(underlying, edge)) == EdgeType::ASCENDING) {
	  const GraphVertex &new_vertex = target(edge, graph);
	  if (get(vertex_type, new_vertex) == unstable_type) {
	    if (lower) { // two maxima with no intersection
	      if (compare(get(distance, *lower), get(distance, new_vertex))) {
		higher = new_vertex;
	      } else {
		higher = *lower;
		lower = new_vertex;
	      }
	    } else {
	      lower = new_vertex;
	    }
	  } else { // intersection
	    higher = intersections_up(graph,
				      vertex_type,
				      underlying,
				      edge_type,
				      new_vertex,
				      intersection);
	  }
	}
      }
      return {*lower, higher};
    }

    template <typename Graph,
	      typename EdgeTypeMap,
	      typename VertexInIterator,
	      typename VertexOutIterator>
    void contour_saddles(const Graph &graph,
			 const EdgeTypeMap &edge_type,
			 VertexInIterator intersection,
			 const VertexInIterator intersections_end,
			 VertexOutIterator saddle_out) {
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      for (; intersection != intersections_end; intersection++) {
	const GraphEdge contour_edge =
	  in_edge_by_type(graph,
			  edge_type,
			  *intersection,
			  EdgeType::CONTOUR);
	*(saddle_out++) = source(contour_edge, graph);
      }
    }

    template <typename Graph,
	      typename EdgeTypeMap>
    void clear_saddle(Graph &master,
		      const EdgeTypeMap &edge_type,
		      const typename boost::graph_traits<Graph>::vertex_descriptor &saddle) {
      using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      std::vector<GraphVertex> intersections;
      intersections.reserve(out_degree(saddle, master));
      for (const GraphEdge &edge : boost::make_iterator_range(out_edges(saddle, master))) {
	if (get(edge_type, edge) == EdgeType::CONTOUR) {
	  const GraphVertex &intersection = target(edge, master);
	  const GraphEdge in = in_edge_by_type(master,
					       edge_type,
					       intersection,
					       EdgeType::ASCENDING);
	  const GraphEdge out = out_edge_by_type(master,
						 edge_type,
						 intersection,
						 EdgeType::ASCENDING);
	  const GraphEdge new_edge =
	    add_edge(source(in, master),
		     target(out, master),
		     master).first;
	  put(edge_type, new_edge, EdgeType::ASCENDING);
	  remove_edge(in, master);
	  remove_edge(out, master);
	  intersections.push_back(intersection);
	}
      }
      clear_vertex(saddle, master);
      for (const GraphVertex &intersection : intersections)
	remove_vertex(intersection, master);
      remove_vertex(saddle, master);
    }

    template <typename Graph,
	      typename KeptSequence,
	      typename IntersectedSequence,
	      typename AddAscendingFn,
	      typename VertexTypeMap,
	      typename EdgeTypeMap>
    void build_back(Graph &master,
		    const KeptSequence &kept,
		    const typename boost::graph_traits<Graph>::vertex_descriptor &higher,
		    const IntersectedSequence &intersected_contours,
		    AddAscendingFn add_ascending,
		    VertexTypeMap &vertex_type,
		    EdgeTypeMap &edge_type) {
      using GraphVertex = typename boost::graph_traits<Graph>::vertex_descriptor;
      using GraphEdge = typename boost::graph_traits<Graph>::edge_descriptor;

      // extend the isolated paths to `higher`
      for (const GraphVertex &vertex : kept) {
	GraphVertex previous_intersection = vertex;
	for (const GraphVertex &contour_vertex : intersected_contours) {
	  const GraphVertex &new_intersection = add_vertex(master);
	  put(vertex_type, new_intersection, VertexType::INTERSECTION);
	  const GraphEdge &ascending_edge =
	    add_ascending(previous_intersection, new_intersection, master).first;
	  put(edge_type, ascending_edge, EdgeType::ASCENDING);
	  const GraphEdge &contour_edge =
	    add_edge(contour_vertex, new_intersection, master).first;
	  put(edge_type, contour_edge, EdgeType::CONTOUR);
	  previous_intersection = new_intersection;
	}
	const GraphEdge &ascending_edge =
	  add_ascending(previous_intersection, higher, master).first;
	put(edge_type, ascending_edge, EdgeType::ASCENDING);
      }
    }

    template <typename Graph>
    std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
    add_directed_edge(const typename boost::graph_traits<Graph>::vertex_descriptor &u,
		      const typename boost::graph_traits<Graph>::vertex_descriptor &v,
		      Graph &graph) {
      return add_edge(u, v, graph);
    }

    template <typename Graph>
    std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
    add_opposite_edge(const typename boost::graph_traits<Graph>::vertex_descriptor &u,
		      const typename boost::graph_traits<Graph>::vertex_descriptor &v,
		      Graph &graph) {
      return add_edge(v, u, graph);
    }
  };

  /*!
    \brief Cancel a saddle vertex of the master graph.
    \tparam Master a model of MutableBidirectionalGraph
    \tparam VertexTypeMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and pebble::VertexType
    as value
    \tparam VertexDistanceMap a model of ReadablePropertyMap with
    boost::graph_traits<Master>::vertex_descriptor as key and a scalar with
    value
    \tparam EdgeTypeMap a model of WritablePropertyMap with
    boost::graph_traits<Master>::edge_descriptor as key and pebble::EdgeType as
    value
    \param master the graph to be modified
    \param vertex_type a property map containing the type of every vertex
    \param distance a property map containing the distance of the centroid and
    the equilibrium or intersection point vertices correspond to
    \param edge_type a property map containing the type of every edge
    \param saddle the saddle vertex to be cancelled
  */
  template <typename Master,
	    typename VertexTypeMap,
	    typename VertexDistanceMap,
	    typename EdgeTypeMap>
  void cancel_saddle(Master &master,
		     const VertexTypeMap &vertex_type,
		     const VertexDistanceMap &distance,
		     const EdgeTypeMap &edge_type,
		     const typename boost::graph_traits<Master>::vertex_descriptor &saddle) {
    using namespace internal;
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using MasterEdge = typename boost::graph_traits<Master>::edge_descriptor;
    using ReversedEdge = typename boost::graph_traits<boost::reverse_graph<Master>>::edge_descriptor;
    using Scalar = typename boost::property_traits<VertexDistanceMap>::value_type;
    using sources = typename boost::inv_adjacency_iterator_generator<Master>::type;

    auto begin = sources(in_edges(saddle, master).first, &master);
    const auto end = sources(in_edges(saddle, master).second, &master);
    boost::container::small_vector<MasterVertex, 2> neighbor_saddles;
    for (; begin != end; ++begin) {
      if (get(vertex_type, *begin) == VertexType::INTERSECTION)
	neighbor_saddles.push_back(source(in_edge_by_type(master,
							  edge_type,
							  *begin,
							  EdgeType::CONTOUR),
					  master));
      else
	neighbor_saddles.push_back(*begin);
    }
    const bool keep_original = neighbor_saddles[0] == neighbor_saddles[1];
    std::vector<MasterVertex> ascending_path;
    MasterVertex lower, higher;
    if (keep_original) {
      std::tie(lower, higher) =
	find_unstables(master,
		       vertex_type,
		       distance,
		       boost::typed_identity_property_map<MasterEdge>(),
		       edge_type,
		       saddle,
		       VertexType::MAX,
		       std::less<Scalar>(),
		       std::back_inserter(ascending_path));
    } else {
      const boost::reverse_graph<Master> reversed(master);
      std::tie(lower, higher) =
	find_unstables(reversed,
		       vertex_type,
		       distance,
		       get(boost::edge_underlying, reversed),
		       edge_type,
		       saddle,
		       VertexType::MIN,
		       std::greater<Scalar>(),
		       std::back_inserter(ascending_path));
    }

    std::vector<MasterVertex> intersected_contours;
    contour_saddles(master,
		    edge_type,
		    ascending_path.cbegin(),
		    ascending_path.cend(),
		    std::back_inserter(intersected_contours));

    std::vector<MasterVertex> descending_path;
    if (keep_original) {
      for (const MasterEdge &edge : boost::make_iterator_range(in_edges(saddle, master)))
	if (get(edge_type, edge) == EdgeType::ASCENDING)
	  intersections_down(master,
			     vertex_type,
			     boost::typed_identity_property_map<MasterEdge>(),
			     edge_type,
			     source(edge, master),
			     std::back_inserter(descending_path));
    } else {
      const boost::reverse_graph<Master> reversed(master);
      for (const ReversedEdge &edge : boost::make_iterator_range(in_edges(saddle, reversed)))
	if (get(edge_type, get(boost::edge_underlying, reversed, edge)) == EdgeType::ASCENDING)
	  intersections_down(reversed,
			     vertex_type,
			     get(boost::edge_underlying, reversed),
			     edge_type,
			     source(edge, reversed),
			     std::back_inserter(descending_path));
    }

    // clear intersections along the isolated paths and
    for (const MasterVertex &vertex : ascending_path) {
      clear_vertex(vertex, master);
      remove_vertex(vertex, master);
    }
    for (const MasterVertex &vertex : descending_path) {
      clear_vertex(vertex, master);
      remove_vertex(vertex, master);
    }
    // clear the saddle
    clear_saddle(master, edge_type, saddle);

    // get the last vertices before `lower`
    // these used to be second to last before the `clear_saddle` call
    std::vector<MasterVertex> kept;
    if (keep_original) {
      kept.reserve(in_degree(lower, master));
      for (const MasterEdge &edge : boost::make_iterator_range(in_edges(lower, master)))
	kept.push_back(source(edge, master));
    } else {
      kept.reserve(out_degree(lower, master));
      for (const MasterEdge &edge : boost::make_iterator_range(out_edges(lower, master)))
	kept.push_back(target(edge, master));
    }

    clear_vertex(lower, master);
    remove_vertex(lower, master);

    if (keep_original) {
      build_back(master,
		 kept,
		 higher,
		 intersected_contours,
		 add_directed_edge<Master>,
		 vertex_type,
		 edge_type);
    } else {
      build_back(master,
		 kept,
		 higher,
		 intersected_contours,
		 add_opposite_edge<Master>,
		 vertex_type,
		 edge_type);
    }
  }
};

#endif // PEBBLE_CANCEL_SADDLE_HPP
