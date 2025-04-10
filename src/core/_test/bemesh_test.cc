//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <gtest/gtest.h>

#include <il/Array.h>
#include <il/Array2D.h>

#include "core/be_mesh.h"
#include "elements/rectangle.h"
#include "elements/segment.h"
#include "elements/triangle.h"
#include "elements/point.h"

/* -------------------------------------------------------------------------- */

// 2D mesh tests
TEST(bemesh_point_2d,pt_1){
    int n_elts = 4;
    il::Array2D<double> coor{n_elts + 1, 2, 0.};
    il::Array2D<double> eltC{n_elts, 2, 0.};
    il::Array2D<il::int_t> conn{n_elts, 1};
    double h = 1;
    for (int i = 0; i < n_elts + 1; i++) {
        coor(i, 0) = i * h;
        coor(i, 1) = i * 2* h;
    }
    for (int i = 0; i < n_elts; i++) {
        conn(i, 0) = i;
    }
    bigwham::BEMesh<bigwham::Point<2>> my_mesh(coor, conn);
    my_mesh.ConstructMesh();
    il::Array2D<double> test = my_mesh.collocation_points();
    bool t = true;
    for (int i = 0; i < n_elts; i++) {
        t = t && ( test(i, 0) == coor(i, 0) ) && ( test(i, 1) == coor(i, 1) ) ;
    }
    ASSERT_TRUE(t);
}


TEST(bemesh_seg, seg_0_1) {
  int n_elts = 4;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<double> eltC{n_elts, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 1;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
    eltC(i, 0) = (coor(i, 0) + coor(i + 1, 0)) / 2.;
  }
  bigwham::BEMesh<bigwham::Segment<0>> my_mesh(coor, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  bool t = true;
  for (int i = 0; i < n_elts; i++) {
    t = t && (abs(test(i, 0) - eltC(i, 0)) < 1.e-6);
  }
  ASSERT_TRUE(t);
}
/* -------------------------------------------------------------------------- */

TEST(bemesh_seg, seg_0_2) {
  int n_elts = 6;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<double> eltC{n_elts, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 0.123;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 1) = i * h;
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
    eltC(i, 0) = (coor(i, 0) + coor(i + 1, 0)) / 2.;
    eltC(i, 1) = (coor(i, 1) + coor(i + 1, 1)) / 2.;
  }
  bigwham::BEMesh<bigwham::Segment<0>> my_mesh(coor, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  bool t = true;
  for (int i = 0; i < n_elts; i++) {
    t = t && (abs(test(i, 1) - eltC(i, 1)) < 1.e-6) &&
        (abs(test(i, 0) - eltC(i, 0)) < 1.e-6);
  }
  ASSERT_TRUE(t);
}
/* -------------------------------------------------------------------------- */

TEST(bemesh_seg, seg_1_1) {
  int n_elts = 4;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<double> eltC{n_elts, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 7.23;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
    eltC(i, 0) = (coor(i, 0) + coor(i + 1, 0)) / 2.;
  }
  bigwham::BEMesh<bigwham::Segment<1>> my_mesh(coor, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  bool t = true;
  for (int i = 0; i < n_elts; i++) { // left col points
    t = t && (abs(test(2 * i, 0) - eltC(i, 0) + (h / 2.) / sqrt(2)) < 1.e-6);
  }
  for (int i = 0; i < n_elts; i++) { // right col points
    t = t &&
        (abs(test(2 * i + 1, 0) - eltC(i, 0) - (h / 2.) / sqrt(2)) < 1.e-6);
  }
  ASSERT_TRUE(t && (test.size(0) == n_elts * 2)); //
}
/* -------------------------------------------------------------------------- */
// Triangular Mesh tests
/* -------------------------------------------------------------------------- */

TEST(bemesh_triangle, triangle_0_1) {
  // mesh consisting of 2 triangles on the e_3 plane
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  xyz(3, 1) = 1.;
  xyz(3, 0) = 1.;
  xyz(3, 2) = 0.;
  il::Array2D<il::int_t> conn{2, 3, 0};
  conn(0, 1) = 1;
  conn(0, 2) = 2;
  conn(1, 0) = 2;
  conn(1, 1) = 1;
  conn(1, 2) = 3;
  bigwham::BEMesh<bigwham::Triangle<0>> my_mesh(xyz, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  ASSERT_TRUE(test.size(0) == my_mesh.num_collocation_points() &&
              abs(test(0, 2)) < 1.e-10 && abs(test(0, 1) - 1. / 3.) < 1.e-10 &&
              abs(test(0, 0) - 1. / 3.) < 1.e-10 &&
              abs(test(1, 0) - 2. / 3.) < 1.e-10 &&
              abs(test(1, 1) - 2. / 3.) < 1.e-10 && abs(test(1, 2)) < 1.e-10);
}
/* -------------------------------------------------------------------------- */

TEST(bemesh_triangle, triangle_0_2) {
  // mesh consisting of 2 triangles at 90 deg from one another
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = 1.;
  xyz(2, 1) = 1.;
  xyz(3, 1) = 1.;
  xyz(3, 0) = 0.;
  xyz(3, 2) = 1.; // at 90 degree
  il::Array2D<il::int_t> conn{2, 3, 0};
  conn(0, 1) = 1;
  conn(0, 2) = 2;
  conn(1, 0) = 2;
  conn(1, 1) = 1;
  conn(1, 2) = 3;
  bigwham::BEMesh<bigwham::Triangle<0>> my_mesh(xyz, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  bigwham::Triangle<0> tri_a;
  il::Array2D<double> xv = my_mesh.vertices(0);
  tri_a.SetElement(xv);
  xv = my_mesh.vertices(1);
  bigwham::Triangle<0> tri_b;
  tri_b.SetElement(xv);
  auto ndots = il::dot(tri_a.normal(), tri_b.normal());
  // std::cout << " normal elt 1::" << tri_a.get_normal()[0] << "-"
  //           << tri_a.get_normal()[1] << "-" << tri_a.get_normal()[2] << "\n";
  // std::cout << " normal elt 2::" << tri_b.get_normal()[0] << "-"
  //           << tri_b.get_normal()[1] << "-" << tri_b.get_normal()[2] << "\n";
  auto ndots_1_1 = il::dot(tri_a.normal(), tri_a.tangent1());
  auto ndots_1_2 = il::dot(tri_a.normal(), tri_a.tangent2());
  auto ndots_2_1 = il::dot(tri_b.normal(), tri_b.tangent1());
  auto ndots_2_2 = il::dot(tri_b.normal(), tri_b.tangent2());
  // std::cout << "ndots" << ndots_2_1 << "\n";
  double eps = 1.e-12;
  ASSERT_TRUE(ndots == 0 && ndots_1_1 == 0 && ndots_1_2 == 0 &&
              abs(ndots_2_1) < eps && abs(ndots_2_2) < eps);
}
/* -------------------------------------------------------------------------- */
// Rectangular mesh
/* -------------------------------------------------------------------------- */

TEST(bemesh_rectangle, rect_0_1) {
  // cartesian mesh with _x*n_y rectangles on the same plane
  il::int_t n_x = 3, n_y = 1;
  il::int_t nelts = n_x * n_y;
  double hx = 1.;
  double hy = 2.;
  il::Array2D<double> coor{(n_x + 1) * (n_y + 1), 3, 0.};
  for (il::int_t k2 = 0; k2 < n_y + 1; k2++) {
    for (il::int_t k = 0; k < n_x + 1; k++) {
      coor(k + k2 * (n_x + 1), 0) = hx * k;
      coor(k + k2 * (n_x + 1), 1) = hy * k2;
    }
  }
  il::Array2D<il::int_t> conn{nelts, 4, 0};
  il::int_t c = 0;
  for (il::int_t k2 = 0; k2 < n_y; k2++) {
    for (il::int_t k = 0; k < n_x; k++) {
      conn(c, 0) = k + k2 * n_y;
      conn(c, 1) = k + 1 + k2 * n_y;
      conn(c, 2) = (k + 1) + k2 * n_y + n_x + 1;
      conn(c, 3) = (k) + k2 * n_y + n_x + 1;
      c++;
    }
  }
  // for (il::int_t c = 0; c < nelts; c++) {
  //   std::cout << " elt #" << c << " conn - " << conn(c, 0) << "-" << conn(c, 1)
  //             << "-" << conn(c, 2) << "-" << conn(c, 3) << "\n";
  // }
  bigwham::BEMesh<bigwham::Rectangle<0>> my_mesh(coor, conn);
  my_mesh.ConstructMesh();
  il::Array2D<double> test = my_mesh.collocation_points();
  bigwham::Rectangle<0> rec_a;
  il::Array2D<double> xv = my_mesh.vertices(0);
  rec_a.SetElement(xv);

  auto aux = il::dot(rec_a.tangent1(), rec_a.normal());
  auto aux2 = il::dot(rec_a.tangent1(), rec_a.tangent2());
  // std::cout << rec_a.size() << "\n";
  ASSERT_TRUE(my_mesh.num_elements() == nelts && aux == 0 && aux2 == 0 &&
              (rec_a.size() == 2.0) &&
              (test.size(0) == my_mesh.num_collocation_points()));
}
