/*!
  \file vertex_dividing.hpp
  \brief Determine which vertices are dividing.

  Let us consider a face \f$F\f$, a halfedge \f$H\f$ in its boundary and the
  vertex \f$V\f$ which is the target of \f$H\f$. \f$V\f$ is <em>dividing with
  respect to \f$F\f$</em> if the line through \f$V\f$ and through the point we
  get by projecting the mesh's centroid onto the plane of \f$F\f$ is intersecting
  the edge of \f$F\f$ opposite to \f$V\f$, dividing it in 2.
*/

#ifndef PEBBLE_VERTEX_DIVIDING_HPP
#define PEBBLE_VERTEX_DIVIDING_HPP 1

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/property_map.h>

namespace pebble {
  /*!
    \brief Determine which vertices are dividing.
    \tparam Mesh a model of FaceListGraph
    \tparam FaceBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3
    scalar values as value
    \tparam DividingPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and bool as value
    \param mesh the input mesh
    \param face_barycentric a property map containing the position of the centroid
    projected to the plane of each face in barycentric coordinates
    \param dividing a property map taking true for a halfedge if the vertex at its
    target is dividing with respect to the face of the halfedge
  */
  template <typename Mesh,
	    typename FaceBarycentricPM,
	    typename DividingPM>
  void vertex_dividing(const Mesh &mesh,
		       const FaceBarycentricPM &face_barycentric,
		       const DividingPM &dividing) {
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using Barycentric = typename boost::property_traits<FaceBarycentricPM>::value_type;

    for (const MeshFace &face : faces(mesh)) {
      const Barycentric &barycentric = get(face_barycentric, face);
      const int negative = (barycentric[0] < 0 ? 1 : 0) +
	(barycentric[1] < 0 ? 1 : 0) + (barycentric[2] < 0 ? 1 : 0);
      if (negative == 2) {
	std::size_t coordinate = 0;
	for (const MeshHalfedge halfedge : halfedges_around_face(halfedge(face, mesh), mesh))
	  put(dividing, halfedge, barycentric[coordinate++] > 0);
      }
    }
  }
};

#endif // PEBBLE_VERTEX_DIVIDING_HPP
