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
#include <il/math.h>

#include "elasticity/fullspace_iso_sp3d_segment/bie_elastostatic_sp3d.h"
#include "elements/boundary_element.h"
/* -------------------------------------------------------------------------- */

TEST(SP3D, test_seg_0_dof_dim) {

  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bigwham::Segment<0> source;
  source.SetElement(xy);
  bigwham::ElasticProperties elas(1, 0.3);
  bigwham::BieElastostatic<bigwham::Segment<0>, bigwham::Segment<0>,
                       bigwham::ElasticKernelType::H>
      test(elas, xy.size(1));
  ASSERT_TRUE(test.dof_dimension() == 2);
}
/* -------------------------------------------------------------------------- */

TEST(SP3D, test_seg_0_dim) {
  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bigwham::Segment<0> source;
  source.SetElement(xy);
  bigwham::ElasticProperties elas(1, 0.3);
  bigwham::BieElastostaticSp3d<bigwham::Segment<0>, bigwham::Segment<0>,
                           bigwham::ElasticKernelType::H>
      test(elas, xy.size(1));
  ASSERT_TRUE(test.spatial_dimension() == 2);
}
/* -------------------------------------------------------------------------- */

TEST(SP3D, test_seg_0_self) {

  il::Array2D<double> xy{2, 2, 0.};
  xy(1, 0) = 1.0;
  bigwham::Segment<0> source;
  source.SetElement(xy);
  bigwham::ElasticProperties elas(1, 0.3);
  bigwham::BieElastostaticSp3d<bigwham::Segment<0>, bigwham::Segment<0>,
                           bigwham::ElasticKernelType::H>
      test(elas, xy.size(1));
  il::Array<double> prop{1, 1000.};
  test.set_kernel_properties(prop);
  // std::cout << "test self effect "
  //           << "\n";
  std::vector<double> test_self = test.influence(source, 0, source, 0);
  // for (int i = 0; i < 4; i++) {
  //   std::cout << test_self[i] << "\n";
  // }
  // old way // would have to be removed at some point !
  il::StaticArray2D<double, 2, 2> xys{0.};
  xys(1, 0) = 1.0;

  // correct values
   il::StaticArray2D<double, 2, 2> stnl{0.};
   stnl(0,0)=0.349791;
   stnl(1,1)=0.349791;

  double epsilon = 1.e-6;
  ASSERT_TRUE(abs(stnl(0, 0) - test_self[0]) < epsilon &&
              abs(stnl(1, 0) - test_self[1]) < epsilon &&
              abs(stnl(0, 1) - test_self[2]) < epsilon &&
              abs(stnl(1, 1) - test_self[3]) < epsilon);
}
/* -------------------------------------------------------------------------- */

TEST(SP3D, test_seg_0_1) {
  il::Array2D<double> xy{2, 2, 0.};
  xy(0, 1) = 2.4;
  xy(1, 0) = 3.0;
  bigwham::Segment<0> source;
  source.SetElement(xy);
  bigwham::ElasticProperties elas(1, 0.3);
  // std::cout << "test on inclined elt "
  //           << "\n";
  il::Array2D<double> xy_r{2, 2, 0.};
  xy_r(0, 0) = 1.0;
  xy_r(1, 0) = 5.0;
  xy_r(0, 1) = 1.0;
  xy_r(1, 1) = 0.0;
  bigwham::Segment<0> receiver;
  receiver.SetElement(xy_r);
  bigwham::BieElastostaticSp3d<bigwham::Segment<0>, bigwham::Segment<0>,
                           bigwham::ElasticKernelType::H>
      test(elas, xy.size(1));
  il::Array<double> prop{1, 1000.};
  test.set_kernel_properties(prop);
  std::vector<double> test_self = test.influence(source, 0, receiver, 0);
  //    for (int i=0;i<3;i++){
  //        std::cout << test_self[i] <<"-"  "\n";
  //    }
  // old way ..... // would have to be removed at some point !
  il::StaticArray2D<double, 2, 2> xys{0.}, xys_r{0.};
  for (int i = 0; i < 2; i++) {
    xys(i, 0) = xy(i, 0);
    xys(i, 1) = xy(i, 1);
    xys_r(i, 0) = xy_r(i, 0);
    xys_r(i, 1) = xy_r(i, 1);
  }

    // correct values from previous code
    il::StaticArray2D<double, 2, 2> stnl{0.};
    stnl(0,0)=-0.0767529;
    stnl(0,1)=0.0719822;
    stnl(1,0)= 0.0678421;
    stnl(1,1)=0.0857816;

  double eps = 3.e-5;
  ASSERT_TRUE(abs(stnl(0, 0) - test_self[0]) < eps &&
              abs(stnl(1, 0) - test_self[1]) < eps &&
              abs(stnl(0, 1) - test_self[2]) < eps &&
              abs(stnl(1, 1) - test_self[3]) < eps);
}
/* -------------------------------------------------------------------------- */
