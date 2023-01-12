/*!
  \file edge_saddle.hpp
  \brief Calculate the projection of the centroid to each of the lines of the
  input mesh's edges.
*/

#ifndef PEBBLE_EDGE_SADDLE_HPP
#define PEBBLE_EDGE_SADDLE_HPP 1

#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>

namespace pebble {
  /*!
    \brief Calculate the projection of the centroid to each of the lines of the
    input mesh's edges.
  
    This function calculates the projection in both Descartes and barycentric coordinates.
    \tparam Mesh a model of EdgeListGraph
    \tparam Point a model of Point_3
    \tparam EdgeMinimumPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and CGAL::Point_3 as value
    \tparam EdgeMinimumBarycentricPM a model of WritablePropertyMap with
    boost::graph_traits<Mesh>::edge_descriptor as key and a Container of 2 scalar
    values as value
    \tparam PointMap a model of ReadablePropertyMap with
    boost::graph_traits<Mesh>::vertex_descriptor as key and CGAL::Point_3 as value
    \param mesh the input mesh
    \param centroid the centroid of the mesh
    \param edge_minimum a property map taking the position of the centroid projected to
    the line of each edge in Descartes coordinates
    \param edge_minimum_barycentric a property map taking the position of the centroid projected to
    the line of each edge in barycentric coordinates
    \param point a property map providing the position of each vertex
  */
  template <typename Mesh,
	    typename Point,
	    typename EdgeMinimumPM,
	    typename EdgeMinimumBarycentricPM,
	    typename PointMap>
  void edge_saddle(const Mesh &mesh,
		   const Point &centroid,
		   const EdgeMinimumPM &edge_minimum,
		   const EdgeMinimumBarycentricPM &edge_minimum_barycentric,
		   const PointMap &point) {
    using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Line = CGAL::Line_3<Kernel>;
    using Plane = CGAL::Plane_3<Kernel>;
    using Barycentric = typename boost::property_traits<EdgeMinimumBarycentricPM>::value_type;

    using namespace CGAL::Barycentric_coordinates;

    for (const MeshEdge &edge : edges(mesh)) {
      const MeshHalfedge he = halfedge(edge, mesh);
      const Point &s = get(point, source(he, mesh));
      const Point &t = get(point, target(he, mesh));
      const Plane normal_plane(centroid, s - t);
      const Line edge_line(s, t);
      const Point minimum = get<Point>(intersection(edge_line, normal_plane).get());
      const Plane edge_plane(centroid, s, t);

      put(edge_minimum,
	  edge,
	  minimum);
      Barycentric barycentric;
      segment_coordinates_2(edge_plane.to_2d(s),
			    edge_plane.to_2d(t),
			    edge_plane.to_2d(minimum),
			    std::begin(barycentric));
      put(edge_minimum_barycentric, edge, barycentric);
    }
  }
}

#endif // PEBBLE_EDGE_SADDLE_HPP
