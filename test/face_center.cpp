#include <gtest/gtest.h>
#include <pebble/face_center.hpp>

#include "tetrahedron.hpp"

using Barycentric = std::array<double, 3>;

TEST(FaceCenterTest, TetrahedronBarycentric) {
  Mesh mesh = tetrahedron();
  const auto face_center_pm =
    mesh.add_property_map<Face, Point>("f:center").first;
  const auto face_barycentric_pm =
    mesh.add_property_map<Face, Barycentric>("f:barycentric").first;
  pebble::face_center(mesh,
		      Point(0, 0, 0),
		      face_center_pm,
		      face_barycentric_pm,
		      mesh.points());
  const Face first = *(mesh.faces().begin());
  const Barycentric coordinates = face_barycentric_pm[first];
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_GE(coordinates[i], 0.0);
    EXPECT_LE(coordinates[i], 1.0);
  }
}
