#ifndef VERTEX_DIVIDING_HPP
#define VERTEX_DIVIDING_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/property_map.h>

namespace pebble {

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
}

#endif // VERTEX_DIVIDING_HPP
