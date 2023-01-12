/*!
  \file axes.hpp
  \brief Compute the axes of a mesh's approximating ellipsoid.
*/

#ifndef PEBBLE_AXES_HPP
#define PEBBLE_AXES_HPP 1

namespace pebble {
  /*!
    \brief Compute the axes of a mesh's approximating ellipsoid.
    \tparam Mesh a model of VertexListGraph
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
    \param mesh the input mesh
    \param point a property map providing the position of each vertex
    \return an array of 3 Vector_3 objects representing the 3 axes of the mesh's approximating ellipsoid from longest to shortest
  */
  template <typename Mesh,
	    typename PointMap>
  std::array<typename CGAL::Kernel_traits<
	       typename boost::property_traits<PointMap>::value_type>::Kernel::Vector_3,
	     3>
  axes(const Mesh &mesh, const PointMap &point) {
    using Point = typename boost::property_traits<PointMap>::value_type;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Scalar = typename Kernel::FT;
    using Vector = typename Kernel::Vector_3;
    using Plane = typename Kernel::Plane_3;
    using Line = typename Kernel::Line_3;
    using VertexIterator = typename boost::graph_traits<Mesh>::vertex_iterator;
    using VertexRange = std::pair<VertexIterator, VertexIterator>;

    const VertexRange range = vertices(mesh);
    
    Scalar max_squared_a = 0, max_squared_b = 0, max_squared_c = 0;
    Vector a, b, c;

    for (VertexIterator v1 = range.first; std::next(v1) != range.second; ++v1) {
      const Point point1 = get(point, *v1);
      for (VertexIterator v2 = std::next(v1); v2 != range.second; ++v2) {
	const Point point2 = get(point, *v2);
	const Scalar squared_distance = CGAL::squared_distance(point1, point2);
	if (squared_distance > max_squared_a) {
	  max_squared_a = squared_distance;
	  a = Vector(point1, point2);
	}
      }
    }

    const Plane projection_plane(CGAL::ORIGIN, a);
    for (VertexIterator v1 = range.first; std::next(v1) != range.second; ++v1) {
      const Point point1 = get(point, *v1);
      const Point projected1 = projection_plane.projection(point1);
      for (VertexIterator v2 = std::next(v1); v2 != range.second; ++v2) {
	const Point point2 = get(point, *v2);
	const Point projected2 = projection_plane.projection(point2);
	const Scalar squared_distance = CGAL::squared_distance(projected1, projected2);
	if (squared_distance > max_squared_b) {
	  max_squared_b = squared_distance;
	  b = Vector(projected1, projected2);
	}
      }
    }

    c = CGAL::cross_product(a, b);
    const Line projection_line(CGAL::ORIGIN, c);
    for (VertexIterator v1 = range.first; std::next(v1) != range.second; ++v1) {
      const Point point1 = get(point, *v1);
      const Point projected1 = projection_line.projection(point1);
      for (VertexIterator v2 = std::next(v1); v2 != range.second; ++v2) {
	const Point point2 = get(point, *v2);
	const Point projected2 = projection_line.projection(point2);
	const Scalar squared_distance = CGAL::squared_distance(projected1, projected2);
	if (squared_distance > max_squared_c)
	  max_squared_c = squared_distance;
      }
    }
    
    c /= CGAL::sqrt(c.squared_length());
    c *= CGAL::sqrt(max_squared_c);

    return {a, b, c};
  }
};

#endif // PEBBLE_AXES_HPP
