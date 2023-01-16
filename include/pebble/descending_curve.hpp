/*!
  \file descending_curve.hpp
  \brief Trace stable-saddle isolated ascending curves.
*/

#ifndef PEBBLE_DESCENDING_CURVE_HPP
#define PEBBLE_DESCENDING_CURVE_HPP 1

#include <CGAL/boost/graph/properties.h>
#include <boost/property_map/property_map.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>

namespace pebble {
  namespace internal {
    /*!
      \brief A single iteration of the recursive isolated ascending curve
      tracing implementation.
    */
    template <typename Mesh, typename Point, typename PointMap,
	      typename FaceMinimumPM, typename FaceBarycentricPM, typename CurveIterator>
    typename boost::graph_traits<Mesh>::face_descriptor
    descending_curve(const Mesh &mesh,
		     const Point &centroid,
		     const PointMap &point,
		     const FaceMinimumPM &face_minimum,
		     const FaceBarycentricPM &face_barycentric,
		     const typename boost::graph_traits<Mesh>::halfedge_descriptor &halfedge,
		     const Point &start,
		     CurveIterator curve) {
      using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
      using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
      using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
      using Barycentric = typename boost::property_traits<FaceBarycentricPM>::value_type;
      using Kernel = typename CGAL::Kernel_traits<Point>::type;
      using Segment = CGAL::Segment_3<Kernel>;
      using Line = CGAL::Line_3<Kernel>;
      using Plane = CGAL::Plane_3<Kernel>;

      const MeshFace &current_face = face(halfedge, mesh);
      const Point &minimum = get(face_minimum, current_face);

      const auto outside = [](const Barycentric &coordinate) {
	return boost::range::count_if(coordinate,
				      [](const typename Kernel::FT &value) {
					return value >= 0 && value <= 1;
				      }) != 3;
      };
      if (outside(get(face_barycentric, current_face))) {
	const Point &source_point = get(point, source(halfedge, mesh));
	const Point &target_point = get(point, target(halfedge, mesh));
	const Point &across_point = get(point, target(next(halfedge, mesh), mesh));
	const Plane divider(across_point, minimum, centroid);
	const bool clockwise = divider.has_on_negative_side(source_point)
	  ? true
	  : (divider.has_on_positive_side(target_point)
	     ? false
	     : divider.has_on_negative_side(start));

	const MeshHalfedge new_halfedge = opposite(clockwise
						   ? next(halfedge, mesh)
						   : prev(halfedge, mesh), mesh);
	const Plane integral(minimum, start, centroid);
	const Point new_start = get<Point>(CGAL::intersection(Line(get(point, source(new_halfedge, mesh)),
								   get(point, target(new_halfedge, mesh))),
							      integral).get());

	*(curve++) = ascending_segment<Mesh, Segment>(current_face, Segment(new_start, start));
	return descending_curve(mesh,
				centroid,
				point,
				face_minimum,
				face_barycentric,
				new_halfedge,
				new_start,
				curve);
      } else {
	*(curve++) = ascending_segment<Mesh, Segment>(current_face, Segment(minimum, start));
	return current_face;
      }
    }
  };

  /*!
    \brief Trace stable-saddle isolated ascending curves for a given range of
    saddles.
    \tparam Mesh a model of FaceGraph
    \tparam Point a model of Point_3
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as
    value
    \tparam EdgeMinimumPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and CGAL::Point_3 as value
    \tparam FaceMinimumPM a model of ReadablePropertyMap
    with boost::graph_traits<Mesh>::face_descriptor as key and CGAL::Point_3 as
    value
    \tparam FaceBarycentricPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a Container of 3
    scalar values as value
    \tparam SaddleIterator a model of InputIterator with
    boost::graph_traits<Mesh>::edge_descriptor as value
    \tparam SegmentPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and a std::forward_list
    of ascending_segment objects as value
    \tparam StablePM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and
    boost::graph_traits<Mesh>::face_descriptor as value
    \param mesh the input mesh
    \param centroid the centroid of the mesh
    \param point a property map providing the position of each vertex
    \param edge_minimum  a property map containing the position of the centroid
    projected to the line of each edge in Descartes coordinates
    \param face_minimum a property map containing the position of the centroid
    projected to the plane of each face in Descartes coordinates
    \param face_barycentric a property map containing the position of
    the centroid projected to the plane of each face in barycentric coordinates
    \param saddles_begin the beginning of the range of edges with saddle
    \param saddles_end the end of the range of edges with saddle
    \param segments a property map storing for every halfedge the segments of
    the stable-saddle isolated ascending curve originating from the saddle on
    the halfedge and descending in the direction of the face adjacent to the halfedge
    \param stable a property map storing for every halfedge the stable point
    at the end of the isolated ascending curve originating from the saddle on the
    halfedge and descending in the direction of the face adjacent to the halfedge
  */
  template <typename Mesh,
	    typename Point,
	    typename PointMap,
	    typename EdgeMinimumPM,
	    typename FaceMinimumPM,
	    typename FaceBarycentricPM,
	    typename SaddleIterator,
	    typename SegmentPM,
	    typename StablePM>
  void descending_curves(const Mesh &mesh,
			 const Point &centroid,
			 const PointMap &point,
			 const EdgeMinimumPM &edge_minimum,
			 const FaceMinimumPM &face_minimum,
			 const FaceBarycentricPM &face_barycentric,
			 SaddleIterator saddle,
			 const SaddleIterator saddles_end,
			 const SegmentPM &segments,
			 const StablePM &stable) {
    using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
    using SegmentSequence = typename boost::property_traits<SegmentPM>::value_type;

    for (; saddle != saddles_end; ++saddle) {
      const MeshHalfedge &h = halfedge(*saddle, mesh);
      for (const MeshHalfedge &halfedge : {h, opposite(h, mesh)}) {
	SegmentSequence &curve = get(segments, halfedge);
	const MeshFace destination =
	  internal::descending_curve(mesh,
				     centroid,
				     point,
				     face_minimum,
				     face_barycentric,
				     halfedge,
				     get(edge_minimum, *saddle),
				     std::front_inserter(curve));
	put(stable, halfedge, destination);
      }
    }
  }
};

#endif // PEBBLE_DESCENDING_CURVE_HPP
