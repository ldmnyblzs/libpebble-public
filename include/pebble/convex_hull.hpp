/*!
   \file convex_hull.hpp
   \brief Read points from an STL file and calculate its convex hull.
*/

#ifndef PEBBLE_CONVEX_HULL_HPP
#define PEBBLE_CONVEX_HULL_HPP 1

#include <vector>
#include <array>
#include <fstream>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/convex_hull_3.h>

#include "types.hpp"

/*!
   \namespace pebble Namespace of the libpebble library.
*/
namespace pebble {
  /*!
     \brief Read points from an STL file and calculate its convex hull.
     \tparam Mesh model of MutableFaceGraph with an internal property map with value type CGAL::vertex_point_t
     \param file path of the input file
     \param mesh output mesh
  */
  template <typename Mesh,
	    typename OriginalFaceMap>
  bool convex_hull(const std::string &file,
		   Mesh &hull,
		   OriginalFaceMap &original) {
    using PointMap = typename boost::property_map<Mesh, CGAL::vertex_point_t>::type;
    using Point = typename boost::property_traits<PointMap>::value_type;
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshEdge = typename boost::graph_traits<Mesh>::edge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;

    Mesh mesh;
    const bool read = CGAL::IO::read_polygon_mesh(file, mesh);
    if (!read)
      return false;
    CGAL::convex_hull_3(mesh, hull);

    // coordinates are copied bit-by-bit so this should work
    const PointMap point = get(CGAL::vertex_point, mesh);
    const PointMap point2 = get(CGAL::vertex_point, hull);
    std::set<std::array<Point, 3>> original_faces;
    for (const MeshFace &face : faces(mesh)) {
      const MeshHalfedge he = halfedge(face, mesh);
      std::array<Point, 3> triangle = {get(point, source(he, mesh)),
	get(point, target(he, mesh)),
	get(point, target(next(he, mesh), mesh))};
      std::sort(triangle.begin(), triangle.end());
      original_faces.insert(triangle);
    }
    for (const MeshFace &face : faces(hull)) {
      const MeshHalfedge he = halfedge(face, hull);
      std::array<Point, 3> triangle = {get(point2, source(he, hull)),
	get(point2, target(he, hull)),
	get(point2, target(next(he, hull), hull))};
      std::sort(triangle.begin(), triangle.end());
      put(original, face, original_faces.contains(triangle));
    }

    bool found;
    do {
      found = false;
      for (const MeshEdge &edge : edges(hull)) {
	const MeshHalfedge he = halfedge(edge, hull);
	const MeshHalfedge op = opposite(he, hull);
	const MeshFace f1 = face(he, hull);
	const MeshFace f2 = face(op, hull);
	if (!get(original, f1) && !get(original, f2)) {
	  std::array<Point, 3> o1 = {
	    get(point2, target(he, hull)),
	    get(point2, target(next(he, hull), hull)),
	    get(point2, target(next(op, hull), hull))
	  };
	  std::sort(o1.begin(), o1.end());
	  std::array<Point, 3> o2 = {
	    get(point2, target(op, hull)),
	    get(point2, target(next(op, hull), hull)),
	    get(point2, target(next(he, hull), hull))
	  };
	  std::sort(o2.begin(), o2.end());
	  if (original_faces.contains(o1) &&
	      original_faces.contains(o2)) {
	    put(original, f1, true);
	    put(original, f2, true);
	    found = true;
	  }
	}
      }
    } while (found);

    return true;
  }
}

#endif // PEBBLE_CONVEX_HULL_HPP
