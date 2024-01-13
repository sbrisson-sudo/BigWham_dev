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
#include "elements/point.h"
#include "elasticity/fullspace_iso_sp3d_segment/bie_elastostatic_sp3d.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_0_influence.h"
#include "elasticity/fullspace_iso_2d_segment/bie_elastostatic_segment_1_influence.h"

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

TEST(TwoDP0,test_several_kernels){
    il::Array2D<double> xy{2, 2, 0.};
    xy(1, 0) = 1.0;
    bie::Segment<0> source{};
    source.SetElement(xy);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,bie::ElasticKernelType::H> testH(elas, xy.size(1));

    std::vector<double> test_self = testH.influence(source, 0, source, 0);

    il::Array2D<double> xy_obs{1, 2, 0.};
    xy_obs(0, 0) = 10.0;
    bie::Point<2> obs_pts;
    obs_pts.SetElement(xy_obs);

    bie::BieElastostatic<bie::Segment<0>, bie::Point<2>,bie::ElasticKernelType::T> test_obsT(elas, xy.size(1));
    auto test_disp = test_obsT.influence(source,0,obs_pts,0);
    std::cout << "test disp size ::" << test_disp.size() <<"\n";
    bie::BieElastostatic<bie::Segment<0>, bie::Point<2>,bie::ElasticKernelType::W> test_obsW(elas, xy.size(1));
    auto test_stress = test_obsW.influence(source,0,obs_pts,0);
    std::cout << "test stress size ::" << test_stress.size() <<"\n";

    ASSERT_TRUE((testH.spatial_dimension() == 2) && (test_obsW.spatial_dimension()==test_obsT.spatial_dimension()));
}


TEST(TwoDP0, test_seg_0_self) {

  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bie::Segment<0> source;
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostaticSp3d<bie::Segment<0>, bie::Segment<0>,
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
    // std::cout << test_self[i] << " " << test_self_sp3d[i] << "\n";
  }
  ASSERT_TRUE(tt);
}
/* -------------------------------------------------------------------------- */

TEST(TwoDP0, test_seg_0_1) {
  il::Array2D<double> xy{2, 2, 0.};
  xy(0, 1) = 2.4;
  xy(1, 0) = 3.0;
  bie::Segment<0> source;
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  // std::cout << "test on inclined elt "
  //           << "\n";
  il::Array2D<double> xy_r{2, 2, 0.};
  xy_r(0, 0) = 1.0;
  xy_r(1, 0) = 5.0;
  xy_r(0, 1) = 1.0;
  xy_r(1, 1) = 0.0;
  bie::Segment<0> receiver;
  receiver.SetElement(xy_r);
  bie::BieElastostatic<bie::Segment<0>, bie::Segment<0>,
                       bie::ElasticKernelType::H>
      test(elas, xy.size(1));

  std::vector<double> test_st = test.influence(source, 0, receiver, 0);

  bie::BieElastostaticSp3d<bie::Segment<0>, bie::Segment<0>,
                           bie::ElasticKernelType::H>
      reftest(elas, xy.size(1));
  il::Array<double> prop{1, 10000.};
  reftest.set_kernel_properties(prop);
  std::vector<double> test_st_sp3d = reftest.influence(source, 0, receiver, 0);

  double epsilon = 1.e-6;
  bool tt = true;
  for (int i = 0; i < 2; i++) {
    tt = tt && (abs(test_st[i] - test_st_sp3d[i]) < epsilon);
  }
  if (!tt) {
    std::cout << "Issue in test" << std::endl;
    for (int i = 0; i < 2; i++) {
      std::cout << std::abs(test_st[i] - test_st_sp3d[i]) << std::endl;
    }
  }

  ASSERT_TRUE(tt);
}
/* -------------------------------------------------------------------------- */
