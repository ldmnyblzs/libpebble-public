/*!
  \file intersections.hpp
  \brief Find intersection points of saddle contour lines and isolated ascending curves.
*/

#ifndef PEBBLE_INTERSECTIONS_HPP
#define PEBBLE_INTERSECTIONS_HPP

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Point_3.h>
#include "types.hpp"

namespace pebble {
  namespace internal {
    template <typename Mesh, typename Kernel>
    struct CloserToCentroid {
      using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using FC = FaceContour<MeshHalfedge, CGAL::Point_3<Kernel>, typename Kernel::FT>;
      using HC = HalfedgeContour<MeshHalfedge, CGAL::Point_3<Kernel>>;
      CGAL::Point_3<Kernel> centroid;
      CloserToCentroid(const CGAL::Point_3<Kernel> &centroid) : centroid(centroid) {}

      bool operator()(const CGAL::Point_3<Kernel> &point,
		      const FC &contour) const {
	return CGAL::has_smaller_distance_to_point(centroid, point, contour.from);
      }
      bool operator()(const FC &contour,
		      const CGAL::Point_3<Kernel> &point) const {
	return CGAL::has_smaller_distance_to_point(centroid, contour.from, point);
      }
      bool operator()(const HC &contour,
		      const CGAL::Point_3<Kernel> &point) const {
	return CGAL::has_smaller_distance_to_point(centroid, contour.at, point);
      }
    };

    template <typename Mesh,
	      typename FaceMinimumPM,
	      typename FaceNormalPM,
	      typename FaceContoursPM,
	      typename Kernel,
	      typename IntersectionIterator>
    void intersections(const Mesh &mesh,
		       const CGAL::Point_3<Kernel> &centroid,
		       const FaceMinimumPM &face_minimum,
		       const FaceNormalPM &face_normal,
		       const FaceContoursPM &face_contours,
		       const typename boost::graph_traits<Mesh>::edge_descriptor &origin,
		       const typename boost::graph_traits<Mesh>::face_descriptor &face,
		       const CGAL::Segment_3<Kernel> &integral,
		       IntersectionIterator intersection) {
      using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using ContourSequence = typename boost::property_traits<FaceContoursPM>::value_type;
      using ContourIterator = typename ContourSequence::const_iterator;
      using Scalar = typename Kernel::FT;
      using Point = CGAL::Point_3<Kernel>;
      using Vector = CGAL::Vector_3<Kernel>;

      const Point &center = get(face_minimum, face);
      const Vector &normal = get(face_normal, face);

      const ContourSequence &contours = get(face_contours, face);
      const ContourIterator lower = std::lower_bound(contours.cbegin(),
						     contours.cend(),
						     integral[0],
						     CloserToCentroid<Mesh, Kernel>(centroid));
      const ContourIterator upper = std::upper_bound(contours.cbegin(),
						     contours.cend(),
						     integral[1],
						     CloserToCentroid<Mesh, Kernel>(centroid));

      for (ContourIterator it = contours.cbegin(); it != lower; it++) {
	if (it->angle == 360 && edge(it->origin, mesh) != origin) {
	  const Scalar integral_angle = angle_on_plane(center, normal, it->from, integral[1]);
	  *(intersection++) =
	    Intersection<MeshHalfedge, Point>{it->origin,
	    center + rotate(it->from - center,
			    normal,
			    integral_angle)};
	}
      }
 
      if (lower < upper) {
	for (ContourIterator it = lower; it != upper; ++it) {
	  if (edge(it->origin, mesh) != origin) {
	    const Scalar integral_angle = angle_on_plane(center, normal, it->from, integral[1]);
	    if (integral_angle <= it->angle + 0.2) {
	      *(intersection++) =
		Intersection<MeshHalfedge, Point>{it->origin,
		center + rotate(it->from - center,
				normal,
				integral_angle)};
	    }
	  }
	}
      }
    }

    template <typename Mesh,
	      typename Kernel,
	      typename HalfedgeContoursPM,
	      typename IntersectionIterator>
    void intersections(const Mesh &mesh,
		       const CGAL::Point_3<Kernel> &centroid,
		       const HalfedgeContoursPM &halfedge_contours,
		       const typename boost::graph_traits<Mesh>::edge_descriptor &origin,
		       const typename boost::graph_traits<Mesh>::halfedge_descriptor &halfedge,
		       const CGAL::Point_3<Kernel> &from,
		       IntersectionIterator intersection) {
      using Point = CGAL::Point_3<Kernel>;
      using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using ContourSequence = typename boost::property_traits<HalfedgeContoursPM>::value_type;
      using ContourIterator = typename ContourSequence::const_iterator;

      const ContourSequence &contours = get(halfedge_contours, halfedge);
      const ContourIterator lower = std::lower_bound(contours.cbegin(),
						     contours.cend(),
						     from,
						     CloserToCentroid<Mesh, Kernel>(centroid));

      if (lower != contours.end()) {
	for (ContourIterator it = lower; it != contours.end(); ++it) {
	  if (edge(it->origin, mesh) != origin) {
	    *(intersection++) = Intersection<MeshHalfedge, Point>{it->origin, it->at};
	  }
	}
      }
    }
  };

  /*!
    \brief Find intersection points of saddle contour lines and isolated
    ascending curves for a given range of saddles.
    \tparam Mesh a model of FaceGraph
    \tparam Point a model of Point_3
    \tparam FaceMinimumPM a model of ReadablePropertyMap
    with boost::graph_traits<Mesh>::face_descriptor as key and CGAL::Point_3 as
    value
    \tparam FaceNormalPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and CGAL::Vector_3 as
    value
    \tparam FaceContoursPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::face_descriptor as key and a std::vector of
    FaceContour objects as value
    \tparam HalfedgeContoursPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and a std::vector of
    HalfedgeContour objects as value
    \tparam SegmentPM a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and a std::vector
    of ascending_segment objects as value
    \tparam SaddleIterator a model of InputIterator with
    boost::graph_traits<Mesh>::edge_descriptor as value
    \tparam IntersectionMap a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::halfedge_descriptor as key and a vector of Intersection
    objects as value
    \param mesh the input mesh
    \param centroid the centroid of the
    mesh
    \param face_minimum a property map containing the position of the centroid
    projected to the plane of each face in Descartes coordinates
    \param face_normal a property map containing the normal vector of each face
    \param face_contours a property map storing for every face the contours it
    was crossed by
    \param halfedge_contours a property map storing for every halfedge the
    contours it was crossed by
    \param segments a property map storing for every halfedge the segments of
    isolated ascending curve originating from the saddle on the halfedge and
    ascending in the direction of the halfedge's target or descending in the
    direction of the face adjacent to the halfedge
    \param saddles_begin the beginning of the range of edges with saddle
    \param saddles_end the end of the range of edges with saddle
    \param intersection_map a property map storing for every halfedge the
    intersection points on the isolated ascending curve originating from the saddle
    on the halfedge
   */
  template <typename Mesh,
	    typename Point,
	    typename FaceMinimumPM,
	    typename FaceNormalPM,
	    typename FaceContoursPM,
	    typename HalfedgeContoursPM,
	    typename SegmentMap,
	    typename SaddleIterator,
	    typename IntersectionMap>
  void intersections(const Mesh &mesh,
		     const Point &centroid,
		     const FaceMinimumPM &face_minimum,
		     const FaceNormalPM &face_normal,
		     const FaceContoursPM &face_contours,
		     const HalfedgeContoursPM &halfedge_contours,
		     const SegmentMap &segment_map,
		     SaddleIterator saddle,
		     const SaddleIterator saddles_end,
		     const IntersectionMap &intersection_map) {
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using AscendingSequence = typename boost::property_traits<SegmentMap>::value_type;
    using AscendingSegment = typename AscendingSequence::value_type;

    for (; saddle != saddles_end; ++saddle) {
      const MeshHalfedge &h = halfedge(*saddle, mesh);
      for (const MeshHalfedge &halfedge : {h, opposite(h, mesh)}) {
	for (const AscendingSegment &segment : get(segment_map, halfedge)) {
	  if (segment.on_face()) {
	    internal::intersections(mesh,
				    centroid,
				    face_minimum,
				    face_normal,
				    face_contours,
				    *saddle,
				    segment.face(),
				    segment.segment(),
				    std::back_inserter(get(intersection_map, halfedge)));
	  } else {
	    internal::intersections(mesh,
				    centroid,
				    halfedge_contours,
				    *saddle,
				    segment.halfedge(),
				    segment.segment()[0],
				    std::back_inserter(get(intersection_map, halfedge)));
	  }
	}
      }
    }
  }
};

#endif // PEBBLE_INTERSECTIONS_HPP
