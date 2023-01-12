/*!
  \file face_center.hpp
  \brief Calculate the projection of the centroid to each of the planes of the input mesh's faces.
 */

#ifndef PEBBLE_FACE_CENTER_HPP
#define PEBBLE_FACE_CENTER_HPP 1

#include <CGAL/boost/graph/properties.h>
#include <boost/property_map/property_map.hpp>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

namespace pebble {
  /*!
    \brief Calculate the projection of the centroid to each of the planes of the
    input mesh's faces.
  
    This function calculates the projection in both Descartes and barycentric coordinates.
    \tparam Mesh a model of FaceListGraph
    \tparam Point a model of Point_3
    \tparam FaceCenterPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and CGAL::Point_3 as value
    \tparam FaceCenterBarycentricPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3 scalar
    values as value
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
    \param mesh the input mesh
    \param centroid the centroid of the mesh
    \param center a property map taking the position of the centroid projected to
    the plane of each face in Descartes coordinates
    \param face_barycentric a property map taking the position of the centroid projected to
    the plane of each face in barycentric coordinates
    \param point a property map providing the position of each vertex
  */
  template <typename Mesh,
	    typename Point,
	    typename FaceCenterPM,
	    typename FaceBarycentricPM,
	    typename PointMap>
  void face_center(const Mesh &mesh,
		   const Point &centroid,
		   const FaceCenterPM &center,
		   const FaceBarycentricPM &face_barycentric,
		   const PointMap &point) {
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Triangle = CGAL::Triangle_3<Kernel>;
    using Plane = CGAL::Plane_3<Kernel>;
    using Scalar = typename Kernel::FT;
    using Barycentric = typename boost::property_traits<FaceBarycentricPM>::value_type;

    using namespace CGAL::Barycentric_coordinates;

    for (const MeshFace &face : faces(mesh)) {
      const MeshHalfedge he = halfedge(face, mesh);
      const Triangle triangle(get(point, source(he, mesh)),
			      get(point, target(he, mesh)),
			      get(point, target(next(he, mesh), mesh)));
      const Plane plane = triangle.supporting_plane();
      const Point projected = plane.projection(centroid);
      put(center, face, projected);
      Barycentric barycentric;
      triangle_coordinates_2(plane.to_2d(triangle[0]),
			     plane.to_2d(triangle[1]),
			     plane.to_2d(triangle[2]),
			     plane.to_2d(projected),
			     std::begin(barycentric));
      put(face_barycentric, face, barycentric);
    }
  }
}

#endif // PEBBLE_FACE_CENTER_HPP
