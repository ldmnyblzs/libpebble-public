#ifndef ASCENDING_CURVE_HPP
#define ASCENDING_CURVE_HPP 1

#include <CGAL/boost/graph/properties.h>
#include <boost/property_map/property_map.hpp>
#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>

namespace pebble {
template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename FaceMinimumPM,
	  typename EdgeMinimumPM,
	  typename EdgeMinimumBarycentricPM,
	  typename IsFollowedPM,
	  typename IsDividingPM,
	  typename CurveIterator>
typename boost::graph_traits<Mesh>::vertex_descriptor
ascending_curve(const Mesh &mesh,
		const Point &centroid,
		const PointMap &point,
		const FaceMinimumPM &face_minimum,
		const EdgeMinimumPM &edge_minimum,
		const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		const IsFollowedPM &is_followed,
		const IsDividingPM &is_dividing,
		const typename boost::graph_traits<Mesh>::vertex_descriptor &vertex,
		CurveIterator curve);

template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename FaceMinimumPM,
	  typename EdgeMinimumPM,
	  typename EdgeMinimumBarycentricPM,
	  typename IsFollowedPM,
	  typename IsDividingPM,
	  typename CurveIterator>
typename boost::graph_traits<Mesh>::vertex_descriptor
ascending_curve(const Mesh &mesh,
		const Point &centroid,
		const PointMap &point,
		const FaceMinimumPM &face_minimum,
		const EdgeMinimumPM &edge_minimum,
		const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		const IsFollowedPM &is_followed,
		const IsDividingPM &is_dividing,
		const typename boost::graph_traits<Mesh>::halfedge_descriptor &halfedge,
		const Point &start,
		CurveIterator curve) {
  using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
  using Barycentric = typename boost::property_traits<EdgeMinimumBarycentricPM>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::type;
  using Segment = CGAL::Segment_3<Kernel>;
  using Line = CGAL::Line_3<Kernel>;
  using Plane = CGAL::Plane_3<Kernel>;

  const MeshFace &current_face = face(halfedge, mesh);
  const MeshEdge &current_edge = edge(halfedge, mesh);
  const Point &minimum = get(face_minimum, current_face);
  const MeshVertex &source_vertex = source(halfedge, mesh);
  const MeshVertex &target_vertex = target(halfedge, mesh);
  const Point &source_point = get(point, source_vertex);
  const Point &target_point = get(point, target_vertex);

  if(get(is_followed, current_edge)) {
    const Barycentric& barycentric = get(edge_minimum_barycentric, current_edge);
    const bool to_target = barycentric[0] > 0 && barycentric[1] > 0
      ? has_smaller_distance_to_point(target_point, start, get(edge_minimum, current_edge))
      : is_ascending(mesh, centroid, halfedge);
    const MeshVertex &vertex = to_target ? target_vertex : source_vertex;
    *(curve++) = ascending_segment<Mesh, Segment>(to_target ? halfedge : opposite(halfedge, mesh),
						  Segment(start, get(point, vertex)));
    return ascending_curve(mesh,
			   centroid,
			   point,
			   face_minimum,
			   edge_minimum,
			   edge_minimum_barycentric,
			   is_followed,
			   is_dividing,
			   vertex,
			   curve);
  } else {
    const Point &across = get(point, target(next(halfedge, mesh), mesh));
    const Plane divider(minimum, across, centroid);

    const bool clockwise = divider.has_on_negative_side(source_point)
      ? true
      : (divider.has_on_positive_side(target_point)
	 ? false
	 : divider.has_on_negative_side(start));
    const MeshHalfedge new_halfedge = opposite(clockwise
					       ? next(halfedge, mesh)
					       : prev(halfedge, mesh), mesh);
    const Plane integral(minimum, start, centroid);
    const Point new_start =
      get<Point>(CGAL::intersection(Line(get(point, source(new_halfedge, mesh)),
					 get(point, target(new_halfedge, mesh))),
				    integral).get());
    *(curve++) = ascending_segment<Mesh, Segment>(current_face, Segment(start, new_start));
    return ascending_curve(mesh,
			   centroid,
			   point,
			   face_minimum,
			   edge_minimum,
			   edge_minimum_barycentric,
			   is_followed,
			   is_dividing,
			   new_halfedge,
			   new_start,
			   curve);
  }
}

template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename FaceMinimumPM,
	  typename EdgeMinimumPM,
	  typename EdgeMinimumBarycentricPM,
	  typename IsFollowedPM,
	  typename IsDividingPM,
	  typename CurveIterator>
typename boost::graph_traits<Mesh>::vertex_descriptor
ascending_curve(const Mesh &mesh,
		const Point &centroid,
		const PointMap &point,
		const FaceMinimumPM &face_minimum,
		const EdgeMinimumPM &edge_minimum,
		const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		const IsFollowedPM &is_followed,
		const IsDividingPM &is_dividing,
		const typename boost::graph_traits<Mesh>::vertex_descriptor &vertex,
		CurveIterator curve) {
  using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
  using Barycentric = typename boost::property_traits<EdgeMinimumBarycentricPM>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::type;
  using Segment = CGAL::Segment_3<Kernel>;
  using Line = CGAL::Line_3<Kernel>;
  using Plane = CGAL::Plane_3<Kernel>;

  const MeshHalfedge he = halfedge(vertex, mesh);
  for (const MeshHalfedge &around : mesh.halfedges_around_target(he)) {
    const MeshHalfedge &opposite_around = opposite(around, mesh);
    const Barycentric &barycentric = get(edge_minimum_barycentric, edge(opposite_around, mesh));
    if (get(is_followed, edge(around, mesh)) &&
	is_ascending(mesh, centroid, opposite_around) &&
	!(barycentric[0] > 0 && barycentric[1] > 0)) {
      *(curve++) =
	ascending_segment<Mesh, Segment>(opposite_around,
					 Segment(get(point, vertex),
						 get(point, target(opposite_around, mesh))));
      return ascending_curve(mesh,
			     centroid,
			     point,
			     face_minimum,
			     edge_minimum,
			     edge_minimum_barycentric,
			     is_followed,
			     is_dividing,
			     source(around, mesh),
			     curve);
    } else if (get(is_dividing, next(around, mesh))) {
      const MeshFace face_around = face(around, mesh);
      const MeshHalfedge across = opposite(next(next(around, mesh), mesh), mesh);
      const Line across_line(get(point, source(across, mesh)), get(point, target(across, mesh)));
      const Plane integral(get(face_minimum, face_around),
			   get(point, target(around, mesh)),
			   centroid);
      const Point new_point = get<Point>(CGAL::intersection(across_line, integral).get());
      *(curve++) = ascending_segment<Mesh, Segment>(face_around,
						    Segment(mesh.point(vertex), new_point));
      return ascending_curve(mesh,
			     centroid,
			     point,
			     face_minimum,
			     edge_minimum,
			     edge_minimum_barycentric,
			     is_followed,
			     is_dividing,
			     across,
			     new_point,
			     curve);
    }
  }
  return vertex;
}

template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename FaceMinimumPM,
	  typename EdgeMinimumPM,
	  typename EdgeMinimumBarycentricPM,
	  typename IsFollowedPM,
	  typename IsDividingPM,
	  typename SaddleIterator,
	  typename SegmentPM,
	  typename UnstablePM>
void ascending_curves(const Mesh &mesh,
		      const Point &centroid,
		      const PointMap &point,
		      const FaceMinimumPM &face_minimum,
		      const EdgeMinimumPM &edge_minimum,
		      const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		      const IsFollowedPM &is_followed,
		      const IsDividingPM &is_dividing,
		      const SaddleIterator &saddles_begin,
		      const SaddleIterator &saddles_end,
		      const SegmentPM &segments,
		      const UnstablePM &unstable) {
  using MeshVertex = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using SegmentSequence = typename boost::property_traits<SegmentPM>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::type;
  using Segment = CGAL::Segment_3<Kernel>;

  for (SaddleIterator saddle = saddles_begin;
       saddle != saddles_end;
       ++saddle) {
    const MeshHalfedge &he = halfedge(*saddle, mesh);
    for (const MeshHalfedge &halfedge : {he, opposite(he, mesh)}) {
      SegmentSequence &curve = get(segments, halfedge);
      curve.emplace_back(halfedge,
			 Segment(get(edge_minimum, *saddle),
				 get(point, target(halfedge, mesh))));
      const MeshVertex &destination =
	ascending_curve(mesh,
			centroid,
			point,
			face_minimum,
			edge_minimum,
			edge_minimum_barycentric,
			is_followed,
			is_dividing,
			target(halfedge, mesh),
			std::back_inserter(curve));
      put(unstable, halfedge, destination);
    }
  }
}

template <typename Mesh,
	  typename Point,
	  typename FaceMinimumPM,
	  typename EdgeMinimumPM,
	  typename EdgeMinimumBarycentricPM,
	  typename IsFollowedPM,
	  typename IsDividingPM,
	  typename SaddleIterator,
	  typename SegmentPM,
	  typename UnstablePM>
void ascending_curves(const Mesh &mesh,
		      const Point &centroid,
		      const FaceMinimumPM &face_minimum,
		      const EdgeMinimumPM &edge_minimum,
		      const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		      const IsFollowedPM &is_followed,
		      const IsDividingPM &is_dividing,
		      const SaddleIterator &saddles_begin,
		      const SaddleIterator &saddles_end,
		      const SegmentPM &segments,
		      const UnstablePM &unstable) {
  ascending_curves(mesh,
		   centroid,
		   get(boost::vertex_point, mesh),
		   face_minimum,
		   edge_minimum,
		   edge_minimum_barycentric,
		   is_followed,
		   is_dividing,
		   saddles_begin,
		   saddles_end,
		   segments,
		   unstable);
}
}

#endif // ASCENDING_CURVE_HPP
