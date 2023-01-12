/*!
  \file types.hpp
  \brief The types storing the mesh and the master graph in pebble-cli.
*/

#ifndef TYPES_HPP
#define TYPES_HPP 1

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <boost/graph/adjacency_list.hpp>
#include <pebble/types.hpp>

using Kernel = CGAL::Simple_cartesian<double>;
using Scalar = Kernel::FT;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Line = Kernel::Line_3;
using Segment = Kernel::Segment_3;
using Plane = Kernel::Plane_3;
using Sphere = Kernel::Sphere_3;
using Mesh = CGAL::Surface_mesh<Point>;
using VertexHandle = Mesh::Vertex_index;
using HalfedgeHandle = Mesh::Halfedge_index;
using EdgeHandle = Mesh::Edge_index;
using FaceHandle = Mesh::Face_index;
using Barycentric = std::array<Scalar, 3>;

using Simplex = std::variant<std::monostate, VertexHandle, EdgeHandle, FaceHandle>;

struct VertexData {
  unsigned int id = 0;
  Point point;
  pebble::VertexType type;
  Simplex simplex;
  Scalar distance;
};

struct EdgeData {
  unsigned int id = 0;
  pebble::EdgeType type;
  bool contour;
};

using Master =
    boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
                          VertexData, EdgeData>;
using MasterVertex = boost::graph_traits<Master>::vertex_descriptor;
using MasterEdge = boost::graph_traits<Master>::edge_descriptor;

#endif // TYPES_HPP
