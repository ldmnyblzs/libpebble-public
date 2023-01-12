#ifndef EDGE_TYPE
#define EDGE_TYPE

#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace pebble {

template <typename Mesh,
	  typename FaceBarycentricPM,
	  typename IsFollowedPM,
	  typename IsNegativePM>
void edge_type(const Mesh &mesh,
	       const FaceBarycentricPM &face_barycentric,
	       const IsFollowedPM &is_followed,
	       IsNegativePM &is_negative) {
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
  using Barycentric = typename boost::property_traits<FaceBarycentricPM>::value_type;
  using HalfedgeIterator = CGAL::Halfedge_around_face_iterator<Mesh>;

  for (const MeshFace &face : faces(mesh)) {
    const Barycentric &barycentric = get(face_barycentric, face);
    HalfedgeIterator he = halfedges_around_face(halfedge(face, mesh), mesh).begin();
    for (const int &index : {2, 0, 1})
      put(is_negative, *(he++), barycentric[index] < 0);
  }

  for (const MeshEdge &edge : edges(mesh)) {
    const MeshHalfedge he = halfedge(edge, mesh);
    put(is_followed,
	edge,
	get(is_negative, he) == get(is_negative, opposite(he, mesh)));
  }

}

template <typename Mesh,
	  typename FaceBarycentricPM,
	  typename IsFollowedPM>
void edge_type(const Mesh &mesh,
	       const FaceBarycentricPM &face_barycentric,
	       const IsFollowedPM &is_followed) {
  using IndexMap = typename boost::property_map<Mesh, CGAL::halfedge_index_t>::type;
  using IsNegativeMap = boost::iterator_property_map<std::vector<bool>::iterator, IndexMap>;

  const IndexMap index = get(CGAL::halfedge_index, mesh);
  std::vector<bool> storage(num_halfedges(mesh));
  IsNegativeMap is_negative_map(storage.begin(), index);
  edge_type(mesh, face_barycentric, is_followed, is_negative_map);
}
}

#endif // EDGE_TYPE
