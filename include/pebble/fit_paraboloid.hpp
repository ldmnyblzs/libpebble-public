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
	    typename FaceOriginalMap,
	    typename Point,
	    typename Master,
	    typename MasterPointMap,
	    typename VertexTypeMap,
	    typename VertexSimplexMap,
	    typename ParaTypeMap>
  void fit_paraboloid(const Mesh &mesh,
		      const MeshPointMap &mesh_point,
		      const FaceOriginalMap &original,
		      const Point &centroid,
		      const Master &master,
		      const MasterPointMap &master_point,
		      const VertexTypeMap &vertex_type,
		      const VertexSimplexMap &simplex,
		      const ParaTypeMap &para_type) {
    using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
    using MasterVertex = typename boost::graph_traits<Master>::vertex_descriptor;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Scalar = typename Kernel::FT;
    using Vector = CGAL::Vector_3<Kernel>;
    using Plane = CGAL::Plane_3<Kernel>;
    using Segment = CGAL::Segment_3<Kernel>;
    using Transformation = CGAL::Aff_transformation_3<Kernel>;
    using SphericalCoordinateSystem =
      boost::geometry::cs::spherical_equatorial<boost::geometry::radian>;
    using SphericalPoint = boost::geometry::model::point<Scalar, 3, SphericalCoordinateSystem>;
    using Matrix3f = Eigen::Matrix<Scalar, 3, 3>;
    using Vector3f = Eigen::Matrix<Scalar, 3, 1>;
    using Matrix6f = Eigen::Matrix<Scalar, 6, 6>;
    using Vector6f = Eigen::Matrix<Scalar, 6, 1>;

    const Transformation translation(CGAL::TRANSLATION, Vector(centroid, Point(0, 0, 0)));

    static const Vector new_origin(1, 0, 0);

    for (const MasterVertex &equilibrium : boost::make_iterator_range(vertices(master)))
      put(para_type, equilibrium, VertexType::INTERSECTION);
    
    std::set<MeshFace> saddle_faces;
    for (const MasterVertex &equilibrium : boost::make_iterator_range(vertices(master))) {
      const Point eq_point = get(master_point, equilibrium);
      if (get(vertex_type, equilibrium) == VertexType::SADDLE) {
	const MeshEdge edge = get<MeshEdge>(get(simplex, equilibrium));
	const MeshHalfedge h1 = halfedge(edge, mesh);
	const MeshHalfedge h2 = opposite(h1, mesh);
	if (!(original[face(h1, mesh)] || original[face(h2, mesh)])) {
	  std::vector<Point> neighborhood;
	  Plane query(eq_point, get(mesh_point, source(h1, mesh))
		      - get(mesh_point, target(h1, mesh)));
	  for (const MeshHalfedge &h : {h1, h2}) {
	    for (const MeshHalfedge &he : CGAL::halfedges_around_face(h, mesh)) {
	      if (he != h) {
		const Segment edge(get(mesh_point, source(he, mesh)),
				 get(mesh_point, target(he, mesh)));
		const auto result = CGAL::intersection(edge, query);
		if (result)
		  if (const Point *point = get<Point>(&*result))
		    neighborhood.push_back(*point);
	      }
	    }
	  }
	  if (neighborhood.size() == 2) {
	    std::array<bool, 2> ascending =
	      {CGAL::has_larger_distance_to_point(centroid, eq_point, neighborhood.at(0)),
	       CGAL::has_larger_distance_to_point(centroid, eq_point, neighborhood.at(1))};
	    if (ascending[0] && ascending[1]) {
	      put(para_type, equilibrium, VertexType::SADDLE);
	      saddle_faces.insert(face(h1, mesh));
	      saddle_faces.insert(face(h2, mesh));
	    } else if (ascending[0] && !ascending[1]) {
	      saddle_faces.insert(face(h2, mesh));
	    } else if (!ascending[0] && ascending[1]) {
	      saddle_faces.insert(face(h1, mesh));
	    } else {
	      const Scalar x1 = -CGAL::approximate_angle(eq_point, centroid, neighborhood.at(0));
	      const Scalar x2 = 0.0;
	      const Scalar x3 = CGAL::approximate_angle(eq_point, centroid, neighborhood.at(1));
	      const Scalar y1 = CGAL::squared_distance(centroid, neighborhood.at(0));
	      const Scalar y2 = CGAL::squared_distance(centroid, eq_point);
	      const Scalar y3 = CGAL::squared_distance(centroid, neighborhood.at(1));
	      const Vector3f b(y1, y2, y3);
	      Matrix3f Am;
	      Am << x1*x1, x1, 1, x2*x2, x2, 1, x3*x3, x3, 1;
	      const Vector3f x = Am.colPivHouseholderQr().solve(b);
	      const Scalar &A = x(0);
	      const Scalar &B = x(1);
	      const Scalar &C = x(2);
	      const Scalar x_min = -B / (2.0 * A);
	      if (x_min > 0.0) {
		saddle_faces.insert(face(h1, mesh));
	      } else {
		saddle_faces.insert(face(h2, mesh));
	      }
	    }
	  }
	}
      }
    }
    
    for (const MasterVertex &equilibrium : boost::make_iterator_range(vertices(master))) {
      const Point eq_point = get(master_point, equilibrium);
      if (get(vertex_type, equilibrium) == VertexType::MIN) {
	const MeshFace face = get<MeshFace>(get(simplex, equilibrium));
	if (!get(original, face) && !saddle_faces.contains(face)) {
	  put(para_type, equilibrium, VertexType::MIN);
	}
      }
    }

    std::map<MasterVertex, std::array<MeshVertex, 6>> unstable_neighborhood;
    for (const MasterVertex &equilibrium : boost::make_iterator_range(vertices(master))) {
      const Point eq_point = get(master_point, equilibrium);
      switch (get(vertex_type, equilibrium)) {
      case VertexType::MIN: {
	const MeshFace face = get<MeshFace>(get(simplex, equilibrium));
	if (!get(original, face))
       	  continue;
       	break;
      }
      case VertexType::SADDLE: {
	const MeshEdge edge = get<MeshEdge>(get(simplex, equilibrium));
	const MeshHalfedge h1 = halfedge(edge, mesh);
	const MeshHalfedge h2 = opposite(h1, mesh);
	if (!(original[face(h1, mesh)] && original[face(h2, mesh)]))
	  continue;
	break;
      }
      case VertexType::MAX: {
	bool skip = true;
	const MeshVertex vertex = get<MeshVertex>(get(simplex, equilibrium));
	for (const MeshFace &face : faces_around_target(halfedge(vertex, mesh), mesh)) {
	  if (original[face]) {
	    skip = false;
	    break;
	  }
	}
	if (skip)
	  continue;
	break;
      }
      }
      if (get(vertex_type, equilibrium) != VertexType::INTERSECTION) {
	// TODO: optimize the algebra
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
	
	std::multimap<Scalar, std::pair<MeshVertex, SphericalPoint>> distances;
	Scalar max_radius = std::numeric_limits<Scalar>::lowest();
	for (const MeshVertex &vertex : boost::make_iterator_range(vertices(mesh))) {
	  const Point position = rotation(translation(get(mesh_point, vertex)));
	  SphericalPoint on_sphere;
	  transform(position, on_sphere);
	  const auto it = mandatory.find(vertex);
	  if (it != mandatory.end() && *it == vertex) {
	    max_radius = std::max(max_radius, spherical_distance(on_sphere));
	    distances.emplace(0.0,
			      std::make_pair(vertex, on_sphere));
	  } else {
	    distances.emplace(spherical_distance(on_sphere),
			      std::make_pair(vertex, on_sphere));
	  }
	}

	std::vector<std::pair<MeshVertex, SphericalPoint>> closest;
	for (auto it = distances.begin();
	     it != distances.end();
	     ++it) {
	  closest.push_back(it->second);
	  if (closest.size() == 6)
	    max_radius = std::max(max_radius, it->first);
	}

	Vector6f b;
	Matrix6f Am;
	for (std::size_t i = 0; i < 6; ++i) {
	  const SphericalPoint &c = closest[i].second;
	  b(i) = get<RADIUS>(c);
	  Am(i, 0) = square(get<THETA>(c));
	  Am(i, 1) = square(get<PHI>(c));
	  Am(i, 2) = get<THETA>(c);
	  Am(i, 3) = get<PHI>(c);
	  Am(i, 4) = get<THETA>(c) * get<PHI>(c);
	  Am(i, 5) = 1;
	}
	const Vector6f x = Am.colPivHouseholderQr().solve(b);

	const Scalar &A = x(0);
	const Scalar &B = x(1);
	const Scalar &C = x(2);
	const Scalar &D = x(3);
	const Scalar &E = x(4);
	const Scalar &F = x(5);

	const Scalar denom = 4.0 * A * B - square(E);
	const Transformation i_translation = translation.inverse();
	const Transformation i_rotation = rotation.inverse();
	if (denom != 0.0) {
	  const Scalar theta = (D * E - 2.0 * B * C) / denom;
	  const Scalar phi = (C * E - 2.0 * A * D) / denom;
	  const Scalar r =
	    A * square(theta) + B * square(phi) + C * theta + D * phi + E * theta * phi + F;

	  SphericalPoint spherical;
	  set<THETA>(spherical, theta);
	  set<PHI>(spherical, phi);
	  set<RADIUS>(spherical, r);
	  
	  const Scalar min_distance = spherical_distance(spherical);
	  if (min_distance < max_radius) {
	    if (denom < 0) {
	      if (get(vertex_type, equilibrium) == VertexType::SADDLE) {
		put(para_type, equilibrium, VertexType::SADDLE);
	      }
	    } else {
	      if (A > 0 && B > 0) {
		if (get(vertex_type, equilibrium) == VertexType::MIN) {
		  SphericalPoint origin;
		  transform(rotation(translation(eq_point)), origin);
		  std::multimap<Scalar, std::pair<MasterVertex, SphericalPoint>> stables;
		  for (const MasterVertex &stable : boost::make_iterator_range(vertices(master))) {
		    if (stable != equilibrium && get(para_type, stable) == VertexType::MIN) {
		      const Point position = rotation(translation(get(master_point, stable)));
		      SphericalPoint on_sphere;
		      transform(position, on_sphere);
		      stables.emplace(spherical_distance(on_sphere),
				      std::make_pair(stable, on_sphere));
		    }
		  }
		  
		  bool doit = true;
		  for (auto it = stables.cbegin();
		       it != stables.cend() && it->first < max_radius;
		       ++it) {
		    if (get<RADIUS>(it->second.second) < get<RADIUS>(origin)) {
		      doit = false;
		    } else {
		      put(para_type, it->second.first, VertexType::INTERSECTION);
		    }
		  }
		  
                  if (doit) {
		    put(para_type, equilibrium, VertexType::MIN);
		  }
		}
	      } else if (A < 0 && B < 0) {
		if (get(vertex_type, equilibrium) == VertexType::MAX) {
		  bool doit = true;

                  std::array<MeshVertex, 6> neighborhood;
		  neighborhood.at(0) = closest.at(0).first;
		  for (std::size_t i = 1; i < 6; ++i) {
		    if (get<RADIUS>(closest.at(i).second) > get<RADIUS>(closest.at(0).second))
		      doit = false;
		    neighborhood.at(i) = closest.at(i).first;
		  }
		  std::sort(neighborhood.begin(), neighborhood.end());
		  
		  for (const auto &[unstable, other] : unstable_neighborhood) {
		    std::vector<MeshVertex> intersection;
		    std::set_intersection(neighborhood.cbegin(),
					  neighborhood.cend(),
					  other.cbegin(),
					  other.cend(),
					  std::back_inserter(intersection));
		    if (!intersection.empty()) {
		      if (CGAL::has_smaller_distance_to_point(centroid,
							      eq_point,
							      get(master_point, unstable))) {
			doit = false;
		      } else {
			put(para_type, unstable, VertexType::INTERSECTION);
		      }
		    }
		  }
		  unstable_neighborhood.emplace(equilibrium, neighborhood);
		    
		  if (doit) {
		    put(para_type, equilibrium, VertexType::MAX);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
};

#endif // PEBBLE_FIT_PARABOLOID
