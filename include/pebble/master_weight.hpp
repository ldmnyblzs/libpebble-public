#ifndef PEBBLE_MASTER_WEIGHT
#define PEBBLE_MASTER_WEIGHT 1

#include <array>
#include <boost/range/iterator_range_core.hpp>
#include "types.hpp"

namespace pebble {
  template <typename Master,
	    typename VertexTypeMap>
  int count_type(const Master &master,
		 const VertexTypeMap &vertex_type,
		 const VertexType &type) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    // TODO: this is not very optimal
    int result = 0;
    for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
      if (get(vertex_type, vertex) == type)
	++result;
    return result;
  }

  template <typename Master,
	    typename VertexTypeMap,
	    typename KeepMap>
  double master_weight(const Master &master,
		       const VertexTypeMap &vertex_type,
		       const KeepMap &keep) {
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    
    int stable = count_type(master, vertex_type, VertexType::MIN);
    int unstable = count_type(master, vertex_type, VertexType::MAX);
    int saddle = stable + unstable - 2;
  
    std::array<int, 3> type_count{0};
    for (const MasterVertex &vertex : boost::make_iterator_range(vertices(master)))
      if (get(vertex_type, vertex) != pebble::VertexType::INTERSECTION &&
	  get(keep, vertex))
	type_count[static_cast<int>(get(vertex_type, vertex))]++;
    
    return static_cast<double>(type_count[0]) / stable +
      static_cast<double>(type_count[1]) / unstable +
      static_cast<double>(type_count[2]) / saddle;
  }
};

#endif // PEBBLE_MASTER_WEIGHT
