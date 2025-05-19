//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>
//
#include "core/be_mesh.h"
#include "elements/segment.h"
#include "hmat/hierarchical_representation.h"
/* -------------------------------------------------------------------------- */

TEST(H_representation, mesh2d_square_1) {
  // simple mesh with close element - such that nothing can be accepted as
  // far-field
  int n_elts = 10;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 0.123;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
  }
  auto my_mesh = std::make_shared<bigwham::BEMesh<bigwham::Segment<0>>>(coor, conn);
  my_mesh->ConstructMesh();
  bool verbose = true;
  auto hr = bigwham::HRepresentationSquareMatrix(my_mesh, 32, 1.0, verbose);
  ASSERT_TRUE(hr->pattern_.n_FRB == 1 && hr->pattern_.n_LRB == 0);
}
/* -------------------------------------------------------------------------- */

TEST(H_representation, mesh2d_square_2) {
  // simple mesh with close element - such that nothing can be accepted as
  // far-field
  int n_elts = 100;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 0.123;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
  }
  auto my_mesh = std::make_shared<bigwham::BEMesh<bigwham::Segment<0>>>(coor, conn);
  my_mesh->ConstructMesh();
  bool verbose = true;
  auto hr = bigwham::HRepresentationSquareMatrix(my_mesh, 50, 0.0, verbose);
  ASSERT_TRUE(hr->pattern_.n_FRB == 4 && hr->pattern_.n_LRB == 0);
}
/* -------------------------------------------------------------------------- */

TEST(H_representation, mesh2d_square_3) {
  // create a mesh with 2 "clusters" of points "far apart" such that the
  // hiearchical representation should have 2 full-block and 2 low-rank block
  int n_elts = 100;
  il::Array2D<double> coor{n_elts + 2, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 0.1;
  for (int i = 0; i < n_elts / 2 + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = n_elts / 2 + 1; i < n_elts + 2; i++) {
    coor(i, 0) = i * h + 60000.;
  }
  for (int i = 0; i < n_elts / 2; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
  }
  for (int i = n_elts / 2; i < n_elts; i++) {
    conn(i, 0) = i + 1;
    conn(i, 1) = i + 2;
  }
  auto my_mesh = std::make_shared<bigwham::BEMesh<bigwham::Segment<0>>>(coor, conn);
  my_mesh->ConstructMesh();
  bool verbose = true;
  auto hr = bigwham::HRepresentationSquareMatrix(my_mesh, n_elts / 2, 4.0, verbose);
  ASSERT_TRUE(hr->pattern_.n_FRB == 2 && hr->pattern_.n_LRB == 2);
}
/* -------------------------------------------------------------------------- */

TEST(H_representation, mesh2d_rectangle_1) {
  // create 2 mesh  one for source, one for receiver representing 2 "clusters"
  // of points "far apart" such that the hiearchical representation should have
  // 1 low-rank block
  int n_elts = 100;
  il::Array2D<double> coor{n_elts + 1, 2, 0.};
  il::Array2D<il::int_t> conn{n_elts, 2};
  double h = 0.01;
  for (int i = 0; i < n_elts + 1; i++) {
    coor(i, 0) = i * h;
  }
  for (int i = 0; i < n_elts; i++) {
    conn(i, 0) = i;
    conn(i, 1) = i + 1;
  }
  auto source_mesh = std::make_shared<bigwham::BEMesh<bigwham::Segment<0>>>(coor, conn);
  source_mesh->ConstructMesh();
  int n_elts_r = 20;
  il::Array2D<double> coor_r{n_elts_r + 1, 2, 0.};
  il::Array2D<il::int_t> conn_r{n_elts_r, 2};

  for (int i = 0; i < n_elts_r + 1; i++) {
    coor_r(i, 0) = i * h * 2;
    coor_r(i, 1) = 4670.;
  }
  for (int i = 0; i < n_elts_r; i++) {
    conn_r(i, 0) = i;
    conn_r(i, 1) = i + 1;
  }
  auto rec_mesh =
      std::make_shared<bigwham::BEMesh<bigwham::Segment<0>>>(coor_r, conn_r);
  rec_mesh->ConstructMesh();
  bool verbose = true;
  auto hr = bigwham::HRepresentationRectangularMatrix(
      source_mesh, rec_mesh, 60, 14.0, verbose);

  ASSERT_TRUE(hr->pattern_.n_FRB == 0 && hr->pattern_.n_LRB == 1);
}
/* -------------------------------------------------------------------------- */
