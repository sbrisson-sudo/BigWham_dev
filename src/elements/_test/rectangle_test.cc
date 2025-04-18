//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 27.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "elements/rectangle.h"

TEST(rectangle, rectangle_0_1) {
  double hx = 1.;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  ASSERT_TRUE(rec0.spatial_dimension() == 3);
}

TEST(rectangle, rectangle_0_2) {
  double hx = 1.3;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  ASSERT_TRUE(rec0.size() == hx * hy);
}

TEST(rectangle, rectangle_0_3) {
  double hx = 1.3;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  auto xc = rec0.centroid();
  ASSERT_TRUE(xc[0] == hx / 2. && xc[1] == hy / 2. && xc[2] == 0.);
}

TEST(rectangle, rectangle_0_4) {
  double hx = 13;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  auto xn = rec0.normal();
  ASSERT_TRUE(xn[0] == 0. && xn[1] == 0. && xn[2] == 1);
}

TEST(rectangle, rectangle_0_5) {
  double hx = 13;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  auto xn = rec0.tangent1();
  ASSERT_TRUE(xn[0] == 1. && xn[1] == 0. && xn[2] == 0);
}

TEST(rectangle, rectangle_0_6) {
  double hx = 13;
  double hy = 2.;
  il::Array2D<double> xyz{4, 3, 0.0};
  xyz(1, 0) = hx;
  xyz(2, 0) = hx;
  xyz(2, 1) = hy;
  xyz(3, 1) = hy;
  bigwham::Rectangle<0> rec0;
  rec0.SetElement(xyz);
  auto xn = rec0.tangent2();
  ASSERT_TRUE(xn[0] == 0. && xn[1] == 1. && xn[2] == 0);
}
