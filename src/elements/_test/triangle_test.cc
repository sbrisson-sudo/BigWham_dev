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

/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_1) {

  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);

  ASSERT_TRUE(tri0.spatial_dimension() == 3);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_2) {
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  ASSERT_TRUE(tri0.num_collocation_points() == 1);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_3) {
  // check rotation matrix
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto R = tri0.rotation_matrix();
  ASSERT_TRUE(R(0, 0) == 1 && R(1, 1) == 1 && R(2, 2) == 1 && R(0, 1) == 0 &&
              R(0, 2) == 0 & R(1, 0) == 0 && R(1, 2) == 0 && R(2, 0) == 0 &&
              R(2, 1) == 0);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_4) {
  // check centroid
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.centroid();
  ASSERT_TRUE(x[0] == 1. / 3 && x[1] == 1. / 3 && x[2] == 0.);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_5) {
  // check nodes
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.nodes();
  ASSERT_TRUE(abs(x(0, 0) - 1. / 3) < 1.e-12 &&
              abs(x(0, 1) - 1. / 3) < 1.e-12 && abs(x(0, 2)) < 1.e-12);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_6) {
  // check vertices
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.vertices();
  ASSERT_TRUE(x(0, 0) == xyz(0, 0) && x(0, 1) == xyz(0, 1) &&
              x(0, 2) == xyz(0, 2));
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_7) {
  // check normal
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.normal();
  ASSERT_TRUE(x[0] == 0 && x[1] == 0 && x[2] == 1);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_8) {
  // check tangent 1
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.tangent1();
  ASSERT_TRUE(x[0] == 1 && x[1] == 0 && x[2] == 0);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_0_9) {
  // check tangent 2
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<0> tri0;
  tri0.SetElement(xyz);
  auto x = tri0.tangent2();
  ASSERT_TRUE(x[0] == 0 && x[1] == 1 && x[2] == 0);
}
/* -------------------------------------------------------------------------- */


TEST(triangle, triangle_0_10) {
    // check tangent 2
    il::Array2D<double> xyz{3, 3, 2.0};
    xyz(0, 0) = 1; xyz(1, 0) = 3.;
    xyz(2, 1) = 4.;
    bie::Triangle<0> tri0;
    tri0.SetElement(xyz);
    il::Array<double> P1{3},P2{3},P3{3};
    for (il::int_t j=0;j<3;j++){
        P1[j]=xyz(0,j)-xyz(0,j);
        P2[j]=xyz(1,j)-xyz(0,j);
        P3[j]=xyz(2,j)-xyz(0,j);
    }
    auto P1a = tri0.ConvertToLocal(P1);
    auto P2a = tri0.ConvertToLocal(P2);
    auto P3a = tri0.ConvertToLocal(P3);
    std::cout <<" P1 local " << P1a[0] <<" |  " << P1a[1] <<" |  " << P1a[2] <<"\n";
    std::cout <<" P2 local " << P2a[0] <<" |  " << P2a[1] <<" |  " << P2a[2] <<"\n";
    std::cout <<" P3 local " << P3a[0] <<" |  " << P3a[1] <<" |  " << P3a[2] <<"\n";

    ASSERT_TRUE(P2a[0] == 2 && P2a[1] == 0 && P2a[2] == 0);
}
/* -------------------------------------------------------------------------- */


TEST(triangle, triangle_2_0) {
  // check # nodes
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.SetElement(xyz);
  auto x = tri2.nodes();
  ASSERT_TRUE(tri2.num_nodes() == 6);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_2_1) {
  //  check dimension
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.SetElement(xyz);
  auto x = tri2.nodes();
  // std::cout << x(0, 0);
  ASSERT_TRUE(tri2.spatial_dimension() == 3);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_2_2) {
  // check # collocation points
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.SetElement(xyz);
  auto x = tri2.nodes();
  // std::cout << x(0, 0);
  ASSERT_TRUE(tri2.num_collocation_points() == 6);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_2_3) {
  // check rotation matrix
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.SetElement(xyz);
  auto R = tri2.rotation_matrix();
  ASSERT_TRUE(R(0, 0) == 1 && R(1, 1) == 1 && R(2, 2) == 1 && R(0, 1) == 0 &&
              R(0, 2) == 0 & R(1, 0) == 0 && R(1, 2) == 0 && R(2, 0) == 0 &&
              R(2, 1) == 0);
}
/* -------------------------------------------------------------------------- */

TEST(triangle, triangle_2_4) {
  // check collocation points location
  il::Array2D<double> xyz{3, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  bie::Triangle<2> tri2;
  tri2.SetElement(xyz);
  auto xcol = tri2.collocation_points();
  ASSERT_TRUE(xcol.size(0) == 6 && xcol.size(1) == 3);
}
