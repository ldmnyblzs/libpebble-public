/*!
  \file measures.hpp
  \brief Measure surface area or volume of a mesh.
 */
#ifndef PEBBLE_MEASURES_HPP
#define PEBBLE_MEASURES_HPP 1

#include <boost/range/numeric.hpp>

namespace pebble {
  namespace internal {
    /*!
      \brief Construct a triangle from the points at the vertices of a triangular
      face.
      \tparam Mesh a model of FaceGraph
      \tparam PointMap a model of ReadablePropertyMap with boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
      \param mesh the input mesh
      \param point a property map providing the position of each vertex
      \param face a face of \p mesh
      \return a CGAL::Triangle_3 constructed from the triangular \p face
    */
    template <typename Mesh, typename PointMap>
    CGAL::Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel::Triangle_3
    face_triangle(const Mesh &mesh,
		  const PointMap &point,
		  const typename boost::graph_traits<Mesh>::face_descriptor face) {
      using Point = typename boost::property_traits<PointMap>::value_type;
      using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
      using Triangle = typename Kernel::Triangle_3;
      using Halfedge = boost::graph_traits<Mesh>::halfedge_descriptor;
  
      const Halfedge h = halfedge(face, mesh);
      return Triangle(get(point, source(h, mesh)), get(point, target(h, mesh)),
		      get(point, target(next(h, mesh), mesh)));
    }
  };

  /*!
    \brief The surface area of a mesh.
    \tparam Mesh a model of FaceListGraph
    \tparam PointMap a model of ReadablePropertyMap with boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
    \param mesh the input mesh
    \param point a property map providing the position of each vertex
  */
  template <typename Mesh, typename PointMap>
  CGAL::Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel::FT
  surface_area(const Mesh &mesh, const PointMap &point) {
    using Point = typename boost::property_traits<PointMap>::value_type;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Scalar = typename Kernel::FT;
    using Face = boost::graph_traits<Mesh>::face_descriptor;

    return boost::accumulate(faces(mesh),
			     0.0,
			     [&](const Scalar &init, const Face &face) {
			       return init
				 + sqrt(internal::face_triangle(mesh, point, face).squared_area());
			     });
  }

  /*!
    \brief The volume of a mesh.
    \tparam Mesh a model of FaceListGraph
    \tparam PointMap a model of ReadablePropertyMap with boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
    \param mesh the input mesh
    \param point a property map providing the position of each vertex
    \param centroid a point in the interior of the mesh
  */
  template <typename Mesh, typename PointMap, typename Point>
  CGAL::Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel::FT
  volume(const Mesh &mesh, const PointMap &point, const Point &centroid) {
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Scalar = typename Kernel::FT;
    using Halfedge = boost::graph_traits<Mesh>::halfedge_descriptor;
    using Face = boost::graph_traits<Mesh>::face_descriptor;
  
    return boost::accumulate(faces(mesh),
			     0.0,
			     [&](const Scalar &init, const Face &face) {
      const Halfedge h = halfedge(face, mesh);
      return init + abs(volume(centroid, get(point, source(h, mesh)),
			       get(point, target(h, mesh)),
			       get(point, target(next(h, mesh), mesh))));
    });
  }
};

#endif // PEBBLE_MEASURES_HPP
