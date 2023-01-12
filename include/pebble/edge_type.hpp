/*!
  \file edge_type.hpp
  \brief Classify edges as followed or crossed.

  For definitions of the edge types see: https://doi.org/10.48550/arXiv.2106.11626
 */

#ifndef PEBBLE_EDGE_TYPE
#define PEBBLE_EDGE_TYPE 1

#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace pebble {
  /*!
    \brief Classify edges as followed or crossed.
    \tparam Mesh model of FaceListGraph
    \tparam FaceBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3 scalar
    values as value
    \tparam IsFollowedPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and bool as value
    \tparam IsNegativePM a model of ReadWritePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and bool as value
    \param mesh the input mesh
    \param face_barycentric a property map containing the position of the centroid
    projected to the plane of each face in barycentric coordinates
    \param is_followed a property map taking true as value for followed edges and
    false for crossed edges
    \param is_negative a utility property map
  */
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

  /*!
    \brief Classify edges as followed or crossed.

    This function calls its 4 parameter overload with a temporary property map as 4th parameter.
    \tparam Mesh model of FaceListGraph
    \tparam FaceBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3 scalar
    values as value
    \tparam FaceBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3 scalar
    values as value
    \tparam IsFollowedPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and bool as value
    \param mesh the input mesh
    \param face_barycentric a property map containing the position of the centroid
    projected to the plane of each face in barycentric coordinates
    \param is_followed a property map taking true as value for followed edges and
    false for crossed edges
  */
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
};

#endif // PEBBLE_EDGE_TYPE
