/*!
   \file convex_hull.hpp
   \brief Read points from an STL file and calculate its convex hull.
*/

#ifndef PEBBLE_CONVEX_HULL_HPP
#define PEBBLE_CONVEX_HULL_HPP 1

#include <fstream>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/IO/STL.h>
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
  template <typename Mesh>
  bool convex_hull(const std::string &file,
		   Mesh &mesh) {
    using PointMap = typename boost::property_map<Mesh, CGAL::vertex_point_t>::type;

    const bool read = CGAL::IO::read_STL(file, mesh);
    if (!read)
      return false;
    const PointMap points = get(CGAL::vertex_point, mesh);
    CGAL::convex_hull_3(boost::begin(points), boost::end(points), mesh);
    return true;
  }
}

#endif // PEBBLE_CONVEX_HULL_HPP
