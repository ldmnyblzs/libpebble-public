#ifndef CONTOUR_LINE_HPP
#define CONTOUR_LINE_HPP

#include <fstream>
#include <CGAL/iterator.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/squared_distance_3.h>
#include <boost/container/static_vector.hpp>
#include <cstdlib>

namespace pebble {

template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename EdgeMinimumPM,
	  typename EdgeBarycentricPM,
	  typename FaceMinimumPM,
	  typename FaceNormalPM,
	  typename Scalar,
	  typename HalfedgeContoursPM,
	  typename FaceContoursPM>
std::pair<Point, typename boost::graph_traits<Mesh>::halfedge_descriptor>
contour_line(const Mesh &mesh,
	     const Point &centroid,
	     const PointMap &point,
	     const EdgeMinimumPM &edge_minimum,
	     const EdgeBarycentricPM &edge_minimum_barycentric,
	     const FaceMinimumPM &face_minimum,
	     const FaceNormalPM &face_normal,
             const typename boost::graph_traits<Mesh>::halfedge_descriptor &origin,
	     const Scalar &distance,
	     const typename boost::graph_traits<Mesh>::halfedge_descriptor &current,
             const Point &previous_point,
	     const HalfedgeContoursPM &halfedge_contours,
	     const FaceContoursPM &face_contours) {
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = CGAL::Vector_3<Kernel>;

  const MeshEdge current_edge = edge(current, mesh);
  if (current_edge == edge(origin, mesh))
    return {previous_point, current};

  Point minimum;
  if (get(edge_minimum_barycentric, current_edge)[0] > 0
      && get(edge_minimum_barycentric, current_edge)[1] > 0)
    minimum = get(edge_minimum, current_edge);
  else if (is_ascending(mesh, centroid, current))
    return contour_line(mesh,
			centroid,
			point,
			edge_minimum,
			edge_minimum_barycentric,
			face_minimum,
			face_normal,
			origin,
			distance,
			next(current, mesh),
			previous_point,
			halfedge_contours,
			face_contours);
  else
    minimum = get(point, target(current, mesh));

  if (squared_distance(centroid, minimum) >= distance)
    return contour_line(mesh,
			centroid,
			point,
			edge_minimum,
			edge_minimum_barycentric,
			face_minimum,
			face_normal,
			origin,
			distance,
			next(current, mesh),
			previous_point,
			halfedge_contours,
			face_contours);

  const Point base = get(edge_minimum, current_edge);
  const Vector direction = normalize(get(point, source(current, mesh)) - base);
  const Scalar length = sqrt(distance - squared_distance(centroid, base));
  const Point new_point = base + (length * direction);

  get(halfedge_contours, opposite(current, mesh)).emplace_back(origin, new_point);

  const MeshFace &current_face = face(current, mesh);
  const Point &center = get(face_minimum, current_face);
  const Vector &normal = get(face_normal, current_face);
  const Scalar angle = angle_on_plane(center, normal, previous_point, new_point);
  get(face_contours, current_face).emplace_back(origin, previous_point, angle);

  return contour_line(mesh,
		      centroid,
		      point,
		      edge_minimum,
		      edge_minimum_barycentric,
		      face_minimum,
		      face_normal,
		      origin,
		      distance,
		      next(opposite(current, mesh), mesh),
		      new_point,
		      halfedge_contours,
		      face_contours);
}

template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename EdgeMinimumPM,
	  typename EdgeBarycentricPM,
	  typename FaceMinimumPM,
	  typename FaceNormalPM,
	  typename EdgeHandle,
	  typename Scalar,
	  typename HalfedgeContoursPM,
	  typename FaceContoursPM>
void contour_line(const Mesh &mesh,
		  const Point &centroid,
		  const PointMap &point,
		  const EdgeMinimumPM &edge_minimum,
		  const EdgeBarycentricPM &edge_minimum_barycentric,
		  const FaceMinimumPM &face_minimum,
		  const FaceNormalPM &face_normal,
		  const EdgeHandle &edge,
		  const Point &saddle,
		  const Scalar &distance,
		  const HalfedgeContoursPM &halfedge_contours,
		  const FaceContoursPM &face_contours) {
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Vector = CGAL::Vector_3<Kernel>;
  using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;

#ifndef BOOST_ASSERT_IS_VOID
  boost::container::static_vector<MeshHalfedge, 2> lasts;
#endif // BOOST_ASSERT_IS_VOID
  
  const MeshHalfedge &h = halfedge(edge, mesh);
  for (const MeshHalfedge &halfedge : {h, opposite(h, mesh)}) {
    const auto [last_point, last_halfedge] =
      contour_line(mesh,
		   centroid,
		   point,
		   edge_minimum,
		   edge_minimum_barycentric,
		   face_minimum,
		   face_normal,
		   halfedge,
		   distance,
		   next(halfedge, mesh),
		   saddle,
		   halfedge_contours,
		   face_contours);
    const MeshFace f = face(last_halfedge, mesh);
    const Point &center = get(face_minimum, f);
    const Vector &normal = get(face_normal, f);
    const Scalar angle = angle_on_plane(center, normal, last_point, saddle);
    get(face_contours, f).emplace_back(halfedge, last_point, angle);
#ifndef BOOST_ASSERT_IS_VOID
    lasts.push_back(last_halfedge);
#endif // BOOST_ASSERT_IS_VOID
  }
#ifndef BOOST_ASSERT_IS_VOID
  BOOST_ASSERT(lasts[0] != lasts[1]);
#endif // BOOST_ASSERT_IS_VOID
}


template <typename Mesh,
	  typename Point,
	  typename PointMap,
	  typename EdgeMinimumPM,
	  typename EdgeBarycentricPM,
	  typename FaceMinimumPM,
	  typename FaceNormalPM,
	  typename SaddleIterator,
	  typename HalfedgeContoursPM,
	  typename FaceContoursPM>
void contour_lines(const Mesh &mesh,
		   const Point &centroid,
		   const PointMap &point,
		   const EdgeMinimumPM &edge_minimum,
		   const EdgeBarycentricPM &edge_minimum_barycentric,
		   const FaceMinimumPM &face_minimum,
		   const FaceNormalPM &face_normal,
		   const SaddleIterator &saddles_begin,
		   const SaddleIterator &saddles_end,
		   const HalfedgeContoursPM &halfedge_contours,
		   const FaceContoursPM &face_contours) {
  using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Scalar = typename Kernel::FT;

  std::vector<std::tuple<Scalar, MeshEdge, Point>> saddles_with_distance;
  for (SaddleIterator edge = saddles_begin;
       edge != saddles_end;
       ++edge) {
    const Point &saddle = get(edge_minimum, *edge);
    saddles_with_distance.emplace_back(squared_distance(centroid, saddle),
				       *edge,
				       saddle);
  }
  std::sort(saddles_with_distance.begin(),
	    saddles_with_distance.end());

  for (const auto &[distance, edge, saddle] : saddles_with_distance)
    contour_line(mesh,
		 centroid,
		 point,
		 edge_minimum,
		 edge_minimum_barycentric,
		 face_minimum,
		 face_normal,
		 edge,
		 saddle,
		 distance,
		 halfedge_contours,
		 face_contours);
}

template <typename Mesh,
	  typename Point,
	  typename EdgeMinimumPM,
	  typename EdgeBarycentricPM,
          typename FaceMinimumPM,
	  typename FaceNormalPM,
          typename SaddleIterator,
	  typename HalfedgeContoursPM,
          typename FaceContoursPM>
void contour_lines(const Mesh &mesh, const Point &centroid,
                   const EdgeMinimumPM &edge_minimum,
		   const EdgeBarycentricPM &edge_minimum_barycentric,
                   const FaceMinimumPM &face_minimum,
                   const FaceNormalPM &face_normal,
                   const SaddleIterator &saddles_begin,
                   const SaddleIterator &saddles_end,
                   const HalfedgeContoursPM &halfedge_contours,
                   const FaceContoursPM &face_contours) {
  contour_lines(mesh,
		centroid,
		get(boost::vertex_point, mesh),
		edge_minimum,
		edge_minimum_barycentric,
                face_minimum,
                face_normal,
                saddles_begin,
                saddles_end,
                halfedge_contours,
                face_contours);
}
}

#endif // CONTOUR_LINE_HPP
