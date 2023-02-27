#include <boost/geometry/core/coordinate_type.hpp>
#ifndef PEBBLE_FIT_PARABOLOID
#define PEBBLE_FIT_PARABOLOID 1

#include <boost/geometry.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/enum.h>
#include <Eigen/Dense>
#include "types.hpp"

#define THETA 0
#define PHI 1
#define RADIUS 2

#define PEBBLE_ADAPT_POINT_TYPE(Scalar, Point)				\
  template<>								\
  struct boost::geometry::traits::tag<Point> {				\
    typedef point_tag type;						\
  };									\
  template<>								\
  struct boost::geometry::traits::coordinate_type<Point> {		\
    typedef Scalar type;						\
  };									\
  template<>								\
  struct boost::geometry::traits::coordinate_system<Point> {		\
    typedef cs::cartesian type;						\
  };									\
  template<>								\
  struct boost::geometry::traits::dimension<Point> : boost::mpl::int_<3> { \
  };									\
  template<>								\
  struct boost::geometry::traits::access<Point, 0> {			\
    static Scalar get(const Point &p) {					\
      return p.x();							\
    }									\
    static void set(Point &p, const Scalar &value) {			\
      p += Point(value, p.y(), p.z()) - p;				\
    }									\
  };									\
  template<>								\
  struct boost::geometry::traits::access<Point, 1> {			\
    static Scalar get(const Point &p) {					\
      return p.y();							\
    }									\
    static void set(Point &p, const Scalar &value) {			\
      p += Point(p.x(), value, p.z()) - p;				\
    }									\
  };									\
  template<>								\
  struct boost::geometry::traits::access<Point, 2> {			\
    static Scalar get(const Point &p) {					\
      return p.z();							\
    }									\
    static void set(Point &p, const Scalar &value) {			\
      p += Point(p.x(), p.y(), value) - p;				\
    }									\
  };

namespace pebble {
  template <typename Scalar>
  inline Scalar square(const Scalar &value) {
    return value * value;
  }

  template <typename SphericalPoint>
  typename boost::geometry::traits::coordinate_type<SphericalPoint>::type
  spherical_distance(const SphericalPoint &point) {
    return square(get<THETA>(point)) + square(get<PHI>(point));
  }
  
  template <typename Mesh,
	    typename MeshPointMap,
	    typename Point,
	    typename Master,
	    typename MasterPointMap,
	    typename VertexTypeMap,
	    typename VertexSimplexMap,
	    typename KeepMap>
  void fit_paraboloid(const Mesh &mesh,
		      const MeshPointMap &mesh_point,
		      const Point &centroid,
		      const Master &master,
		      const MasterPointMap &master_point,
		      const VertexTypeMap &vertex_type,
		      const VertexSimplexMap &simplex,
		      const KeepMap &keep) {
    using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Scalar = typename Kernel::FT;
    using Vector = CGAL::Vector_3<Kernel>;
    using Transformation = CGAL::Aff_transformation_3<Kernel>;
    using SphericalCoordinateSystem =
      boost::geometry::cs::spherical_equatorial<boost::geometry::radian>;
    using SphericalPoint = boost::geometry::model::point<Scalar, 3, SphericalCoordinateSystem>;
    using Matrix6f = Eigen::Matrix<Scalar, 6, 6>;
    using Vector6f = Eigen::Matrix<Scalar, 6, 1>;

    const Transformation translation(CGAL::TRANSLATION, Vector(centroid, Point(0, 0, 0)));

    static const Vector new_origin(1, 0, 0);
    for (const MasterVertex &equilibrium : boost::make_iterator_range(vertices(master))) {
      put(keep, equilibrium, false);
      if (get(vertex_type, equilibrium) != VertexType::INTERSECTION) {
	// TODO: optimize the algebra
	const Point eq_point = get(master_point, equilibrium);
	const Vector old_origin = eq_point - centroid;
	const Vector axis = normalize(CGAL::cross_product(old_origin, new_origin));
	const Scalar angle = to_radians(CGAL::approximate_angle(old_origin, new_origin));
	const Scalar sinus = sin(angle);
	const Scalar cosinus = cos(angle);
	const Scalar mincos = 1.0 - cosinus;
	const Transformation rotation(cosinus + square(axis.x()) * mincos,
				      axis.x() * axis.y() * mincos - axis.z() * sinus,
				      axis.x() * axis.z() * mincos + axis.y() * sinus,
				      axis.y() * axis.x() * mincos + axis.z() * sinus,
				      cosinus + square(axis.y()) * mincos,
				      axis.y() * axis.z() * mincos - axis.x() * sinus,
				      axis.z() * axis.x() * mincos - axis.y() * sinus,
				      axis.z() * axis.y() * mincos + axis.x() * sinus,
				      cosinus + square(axis.z()) * mincos);

	std::set<MeshVertex> mandatory;
	switch (get(vertex_type, equilibrium)) {
	case VertexType::MIN: {
	  const MeshFace face = get<MeshFace>(get(simplex, equilibrium));
	  for (const MeshVertex &v : vertices_around_face(halfedge(face, mesh), mesh)) {
	    mandatory.insert(v);
	  }
	  break;
	}
	case VertexType::SADDLE: {
	  const MeshEdge edge = get<MeshEdge>(get(simplex, equilibrium));
	  const MeshHalfedge h = halfedge(edge, mesh);
	  for (const MeshVertex v : {source(h, mesh), target(h, mesh)}) {
	    mandatory.insert(v);
	  }
	  break;
	}
	case VertexType::MAX: {
	  mandatory.insert(get<MeshVertex>(get(simplex, equilibrium)));
	  break;
	}
	}
	std::multimap<Scalar, SphericalPoint> distances;
	Scalar max_radius = std::numeric_limits<Scalar>::lowest();
	for (const MeshVertex &vertex : boost::make_iterator_range(vertices(mesh))) {
	  const Point position = rotation(translation(get(mesh_point, vertex)));
	  SphericalPoint on_sphere;
	  transform(position, on_sphere);
	  if (mandatory.find(vertex) != mandatory.end()) {
	    max_radius = std::max(max_radius, spherical_distance(on_sphere));
	    distances.emplace(0.0, on_sphere);
	  } else {
	    distances.emplace(spherical_distance(on_sphere), on_sphere);
	  }
	}

	typename std::multimap<Scalar, SphericalPoint>::const_iterator it =
	  distances.begin();
	std::array<SphericalPoint, 6> closest;
	for (std::size_t i = 0; i < 6; ++i) {
	  max_radius = std::max(max_radius, it->first);
	  closest[i] = it->second;
	  it++;
	}

	Vector6f b;
	Matrix6f Am;
	for (std::size_t i = 0; i < 6; ++i) {
	  b(i) = get<RADIUS>(closest[i]);
	  Am(i, 0) = square(get<THETA>(closest[i]));
	  Am(i, 1) = square(get<PHI>(closest[i]));
	  Am(i, 2) = get<THETA>(closest[i]) * get<PHI>(closest[i]);
	  Am(i, 3) = get<THETA>(closest[i]);
	  Am(i, 4) = get<PHI>(closest[i]);
	  Am(i, 5) = 1;
	}
	const Vector6f x = Am.colPivHouseholderQr().solve(b);

	const Scalar &A = x(0);
	const Scalar &B = x(1);
	const Scalar &C = x(2);
	const Scalar &D = x(3);
	const Scalar &E = x(4);
	const Scalar &F = x(5);

	const Scalar denom = square(C) - 4.0 * A * B;
	if (denom != 0.0 && C != 0.0) {
	  const Scalar phi = (2.0 * A * E - C * D) / denom;
	  const Scalar theta = (-2.0 * B * phi - E) / C;
	  const Scalar r =
	    A * square(theta) + B * square(phi) + C * theta * phi + D * theta + E * phi + F;
	  const Scalar min_distance = square(theta) + square(phi);
	  if (min_distance < max_radius) {
	    const bool Aless = A < 0;
	    const bool Bless = B < 0;
	    switch (get(vertex_type, equilibrium)) {
	    case VertexType::MIN:
	      if (!Aless && !Bless)
		put(keep, equilibrium, true);
	      break;
	    case VertexType::SADDLE:
	      if (Aless != Bless)
		put(keep, equilibrium, true);
	      break;
	    case VertexType::MAX:
	      if (Aless && Bless)
		put(keep, equilibrium, true);
	      break;
	    }
	  }
	}
      }
    }
  }
};

#endif // PEBBLE_FIT_PARABOLOID
