//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 10.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "core/elastic_properties.h"
#include "elasticity/bie_elastostatic.h"
#include "elements/segment.h"
/* -------------------------------------------------------------------------- */

TEST(TwoDP0, test_seg_0_dof_dim) {

  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bie::Segment<0> source{};
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
                       bie::ElasticKernelType::H>
      test(elas, 2);
  ASSERT_TRUE(test.dof_dimension() == 2);
}
/* -------------------------------------------------------------------------- */

TEST(TwoDP0, test_seg_0_dim) {
  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bie::Segment<0> source{};
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
                       bie::ElasticKernelType::H>
      test(elas, xy.size(1));
  ASSERT_TRUE(test.spatial_dimension() == 2);
}
/* -------------------------------------------------------------------------- */

TEST(TwoDP0, test_seg_0_self) {

  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bie::Segment<0> source;
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
                       bie::ElasticKernelType::H>
      reftest(elas, xy.size(1));
  il::Array<double> prop{1, 1000.};
  reftest.set_kernel_properties(prop);
  std::vector<double> test_self_sp3d = reftest.influence(source, 0, source, 0);

  bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
                       bie::ElasticKernelType::H>
      test(elas, xy.size(1));
  std::vector<double> test_self = test.influence(source, 0, source, 0);

  double epsilon = 1.e-6;
  bool tt = true;
  for (int i = 0; i < 3; i++) {
    tt = tt && (abs(test_self[i] - test_self_sp3d[i]) < epsilon);
  }
  ASSERT_TRUE(tt);
}
/* -------------------------------------------------------------------------- */

// TEST(TwoDP0, test_seg_0_1) {
//   il::Array2D<double> xy{2, 2, 0.};
//   xy(0, 1) = 2.4;
//   xy(1, 0) = 3.0;
//   bie::Segment<0> source;
//   source.SetElement(xy);
//   bie::ElasticProperties elas(1, 0.3);
//   std::cout << "test on inclined elt "
//             << "\n";
//   il::Array2D<double> xy_r{2, 2, 0.};
//   xy_r(0, 0) = 1.0;
//   xy_r(1, 0) = 5.0;
//   xy_r(0, 1) = 1.0;
//   xy_r(1, 1) = 0.0;
//   bie::Segment<0> receiver;
//   receiver.SetElement(xy_r);
//   bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
//                        bie::ElasticKernelType::H>
//       test(elas, xy.size(1));

//   std::vector<double> test_st = test.influence(source, 0, receiver, 0);

//   bie::BieElastostaticSp3d<bie::Segment<0>, bie::Segment<0>,
//                            bie::ElasticKernelType::H>
//       reftest(elas, xy.size(1));
//   il::Array<double> prop{1, 10000.};
//   reftest.setKernelProperties(prop);
//   std::vector<double> test_st_sp3d = reftest.influence(source, 0, receiver, 0);

//   double epsilon = 1.e-6;
//   bool tt = true;
//   for (int i = 0; i < 2; i++) {
//     tt = tt && (abs(test_st[i] - test_st_sp3d[i]) < epsilon);
//   }
//   if (!tt) {
//     std::cout << "Issue in test" << std::endl;
//     for (int i = 0; i < 2; i++) {
//       std::cout << std::abs(test_st[i] - test_st_sp3d[i]) << std::endl;
//     }
//   }

//   ASSERT_TRUE(tt);
// }
/* -------------------------------------------------------------------------- */
