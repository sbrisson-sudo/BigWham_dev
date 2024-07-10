//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2024.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <iostream>

#include <gtest/gtest.h>

#include <il/StaticArray2D.h>

#include "elasticity/fullspace_iso_2d_segment/elastic_2dP0_segment.h"

TEST(elas_fs_iso_segment, U_0_1) {

  double h = 2.0;
  double G = 1.5;
  double nu = 0.2;
  double x1 = 0.;
  double x2 = -1.3;
  auto Ue = bigwham::Ue_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 2> U_exact{0.};
  U_exact(0, 0) = -0.0408623;
  U_exact(1, 1) = 0.00587676;
  double eps = 1.e-6;
  // std::cout << Ue(0, 0) << "\n";
  ASSERT_TRUE(abs(Ue(0, 0) - U_exact(0, 0)) < eps &&
              abs(Ue(1, 0) - U_exact(1, 0)) < eps &&
              abs(Ue(1, 1) - U_exact(1, 1)) < eps &&
              abs(Ue(0, 1) - U_exact(0, 1)) < eps);
}

TEST(elas_fs_iso_segment, U_0_2) {
  // test origin self effect
  double h = 2.0;
  double G = 10.5;
  double nu = 0.2;
  double x1 = 0.;
  double x2 = 0.;
  auto Ue = bigwham::Ue_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 2> U_exact{0.};
  U_exact(0, 0) = 0.0303152;
  U_exact(1, 1) = 0.0208417;
  double eps = 1.e-6;
  ASSERT_TRUE(abs(Ue(0, 0) - U_exact(0, 0)) < eps &&
              abs(Ue(1, 0) - U_exact(1, 0)) < eps &&
              abs(Ue(1, 1) - U_exact(1, 1)) < eps &&
              abs(Ue(0, 1) - U_exact(0, 1)) < eps);
}

TEST(elas_fs_iso_segment, U_0_3) {
  //
  double h = 2.0;
  double G = 17.25;
  double nu = 0.22;
  double x1 = 0.2;
  double x2 = 0.8;
  auto Ue = bigwham::Ue_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 2> U_exact{0.};
  // {{0.00202768, 0.000574292}, {0.000574292, 0.0044795}}
  U_exact(0, 0) = 0.00202768;
  U_exact(0, 1) = 0.000574292;
  U_exact(1, 0) = 0.000574292;
  U_exact(1, 1) = 0.0044795;

  double eps = 1.e-6;
  ASSERT_TRUE(abs(Ue(0, 0) - U_exact(0, 0)) < eps &&
              abs(Ue(1, 0) - U_exact(1, 0)) < eps &&
              abs(Ue(1, 1) - U_exact(1, 1)) < eps &&
              abs(Ue(0, 1) - U_exact(0, 1)) < eps);
}

TEST(elas_fs_iso_segment, U_0_4) {
  //
  double h = 2.0;
  double G = 17.25;
  double nu = 0.22;
  double x1 = 0.8;
  double x2 = -0.2;
  auto Ue = bigwham::Ue_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 2> U_exact{0.};
  // {{0.0091904, -0.00109817}, {-0.00109817, 0.00593223}}
  U_exact(0, 0) = 0.0091904;
  U_exact(0, 1) = -0.00109817;
  U_exact(1, 0) = -0.00109817;
  U_exact(1, 1) = 0.00593223;
  double eps = 1.e-6;
  ASSERT_TRUE(abs(Ue(0, 0) - U_exact(0, 0)) < eps &&
              abs(Ue(1, 0) - U_exact(1, 0)) < eps &&
              abs(Ue(1, 1) - U_exact(1, 1)) < eps &&
              abs(Ue(0, 1) - U_exact(0, 1)) < eps);
}

TEST(elas_fs_iso_segment, S_0_1) {
  //
  double h = 2.0;
  double G = 17.25;
  double nu = 0.2;
  double x1 = 0.3;
  double x2 = 8.0;
  auto Se = bigwham::Se_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 3> S_exact{0.};
  //  {{{-0.000580783}, {-0.0151433}}, {{-0.0151433}, {-0.00125313}}}
  S_exact(0, 0) = -0.000580783;
  S_exact(0, 1) = -0.0151433;
  S_exact(0, 2) = -0.00125313;
  //{{{0.0145037}, {-0.00235348}}, {{-0.00235348}, {-0.0639152}}}
  S_exact(1, 0) = 0.0145037;
  S_exact(1, 1) = -0.00235348;
  S_exact(1, 2) = -0.0639152;
  double eps = 1.e-6;
  bool test = true;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      test = test && (abs(Se(i, j) - S_exact(i, j)) < eps);
    }
  }
  ASSERT_TRUE(test);
}

TEST(elas_fs_iso_segment, S_0_2) {
  // on x2=0
  double h = 2.0;
  double G = 17.25;
  double nu = 0.2;
  double x1 = 10.3;
  double x2 = 0;
  auto Se = bigwham::Se_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 3> S_exact{0.};
  //
  S_exact(0, 0) = -0.0503775;
  S_exact(0, 1) = 0.;
  S_exact(0, 2) = 0.0116256;
  S_exact(1, 0) = 0.;
  S_exact(1, 1) = -0.0116256;
  S_exact(1, 2) = 0.;
  double eps = 1.e-3;
  bool test = true;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
        std::cout <<  Se(i, j)  << "\n";
      test = test && (abs(Se(i, j) - S_exact(i, j)) < eps);
    }
  }
  ASSERT_TRUE(test);
}

// Tests for integtation of Vijk2 on the segment [-1,1]
TEST(elas_fs_iso_segment, W_0_1) {
  // on x2=0
  double h = 2.0;
  double G = 17.25;
  double nu = 0.2;
  double x1 = 10.3;
  double x2 = 0;
  auto We = bigwham::We_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 3> W_exact{0.};
  //
  W_exact(0, 0) = 0.;
  W_exact(0, 1) = -0.0653112;
  W_exact(0, 2) = 0.;
  W_exact(1, 0) = -0.0653112;
  W_exact(1, 1) = 0.;
  W_exact(1, 2) = -0.0653112;

  double eps = 1.e-6;
  bool test = true;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      test = test && (abs(We(i, j) - W_exact(i, j)) < eps);
    }
  }
  ASSERT_TRUE(test);
}

// Tests for integtation of Vijk2 on the segment [-1,1]
TEST(elas_fs_iso_segment, W_0_2) {
  // on x2=0
  double h = 2.0;
  double G = 17.25;
  double nu = 0.2;
  double x1 = 1.3;
  double x2 = -2.4;
  auto We = bigwham::We_segment_0(h, G, nu, x1, x2);
  il::StaticArray2D<double, 2, 3> W_exact{0.};
  //
  W_exact(0, 0) = 0.0077946525229993305;
  W_exact(0, 1) = 0.1401913225410614;
  W_exact(0, 2) = -1.3328855814328973;
  W_exact(1, 0) = 0.1401913225410614;
  W_exact(1, 1) = -1.3328855814328973;
  W_exact(1, 2) = 0.9364450571982305;

  double eps = 1.e-11;
  bool test = true;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      test = test && (abs(We(i, j) - W_exact(i, j)) < eps);
    }
  }
  ASSERT_TRUE(test);
}
