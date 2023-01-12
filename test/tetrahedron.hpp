#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Vertex = Mesh::Vertex_index;
using Face = Mesh::Face_index;

Mesh tetrahedron() {
  Mesh mesh;
  const std::array<Vertex, 4> vertices = {
    mesh.add_vertex(Point(1, 1, 1)),
    mesh.add_vertex(Point(-1, 1, -1)),
    mesh.add_vertex(Point(1, -1, -1)),
    mesh.add_vertex(Point(-1, -1, 1))};
  const static std::array<std::array<std::size_t, 3>, 4> faces =
    {{{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}};
  for (const auto &[i, j, k] : faces)
    mesh.add_face(vertices.at(i),
		  vertices.at(j),
		  vertices.at(k));
  return mesh;
}
