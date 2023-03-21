//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "elements/triangle.h"

TEST(triangle, triangle_0_1) {

  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);

  ASSERT_TRUE(tri0.get_spatial_dimension() == 3);
}

TEST(triangle, triangle_0_2) {
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  ASSERT_TRUE(tri0.get_number_collocation_points() == 1);
}

TEST(triangle, triangle_0_3) {
  // check rotation matrix
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto R = tri0.get_rotation_matrix();
  ASSERT_TRUE(R(0, 0) == 1 && R(1, 1) == 1 && R(2, 2) == 1 && R(0, 1) == 0 &&
              R(0, 2) == 0 & R(1, 0) == 0 && R(1, 2) == 0 && R(2, 0) == 0 &&
              R(2, 1) == 0);
}

TEST(triangle, triangle_0_4) {
  // check centroid
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_centroid();
  ASSERT_TRUE(x[0] == 1. / 3 && x[1] == 1. / 3 && x[2] == 0.);
}

TEST(triangle, triangle_0_5) {
  // check nodes
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_nodes();
  ASSERT_TRUE(abs(x(0, 0) - 1. / 3) < 1.e-12 &&
              abs(x(0, 1) - 1. / 3) < 1.e-12 && abs(x(0, 2)) < 1.e-12);
}

TEST(triangle, triangle_0_6) {
  // check vertices
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_vertices();
  ASSERT_TRUE(x(0, 0) == xyz(0, 0) && x(0, 1) == xyz(0, 1) &&
              x(0, 2) == xyz(0, 2));
}

TEST(triangle, triangle_0_7) {
  // check normal
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_normal();
  ASSERT_TRUE(x[0] == 0 && x[1] == 0 && x[2] == 1);
}

TEST(triangle, triangle_0_8) {
  // check tangent 1
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_tangent_1();
  ASSERT_TRUE(x[0] == 1 && x[1] == 0 && x[2] == 0);
}

TEST(triangle, triangle_0_9) {
  // check tangent 2
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.set_element(xyz);
  auto x = tri0.get_tangent_2();
  ASSERT_TRUE(x[0] == 0 && x[1] == 1 && x[2] == 0);
}

TEST(triangle, triangle_2_0) {
  // check # nodes
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.set_element(xyz);
  auto x = tri2.get_nodes();
  ASSERT_TRUE(tri2.get_number_nodes() == 6);
}

TEST(triangle, triangle_2_1) {
  //  check dimension
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.set_element(xyz);
  auto x = tri2.get_nodes();
  std::cout << x(0, 0);
  ASSERT_TRUE(tri2.get_spatial_dimension() == 3);
}

TEST(triangle, triangle_2_2) {
  // check # collocation points
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.set_element(xyz);
  auto x = tri2.get_nodes();
  std::cout << x(0, 0);
  ASSERT_TRUE(tri2.get_number_collocation_points() == 6);
}

TEST(triangle, triangle_2_3) {
  // check rotation matrix
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.set_element(xyz);
  auto R = tri2.get_rotation_matrix();
  ASSERT_TRUE(R(0, 0) == 1 && R(1, 1) == 1 && R(2, 2) == 1 && R(0, 1) == 0 &&
              R(0, 2) == 0 & R(1, 0) == 0 && R(1, 2) == 0 && R(2, 0) == 0 &&
              R(2, 1) == 0);
}

TEST(triangle, triangle_2_4) {
  // check collocation points location
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.set_element(xyz);
  auto xcol = tri2.get_collocation_points();
  ASSERT_TRUE(xcol.size(0) == 6 && xcol.size(1) == 3);
}
