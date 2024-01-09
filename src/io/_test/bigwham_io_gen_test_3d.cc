//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 31.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "hmat/hmatrix/hmat.h"
#include "io/bigwham_io_gen.h"

#include "core/be_mesh.h"
#include "elements/boundary_element.h"
#include "hmat/square_matrix_generator.h"
/* -------------------------------------------------------------------------- */

TEST(bigwham_io_gen_3d, 3DT0_1) {
    // 2 Triangle for the square [0,1]X[0,1]
    int n_elts = 2;
    std::vector<double> coor(3*4  , 0.);
    coor[3]=1.;
    coor[7]=1.; coor[8]=0.;
    coor[9]=1.;coor[10]=1;
    std::vector<int> conn(n_elts * 3, 0);
    int k = 0;
    for (int i = 0; i < n_elts; i++) {
        conn[k] = i;
        conn[k + 1] = i + 1;
        conn[k + 2] = i + 2;
        k = k + 3;
    }
  std::vector<double> properties{1., 0.};

  BigWhamIOGen my_io{coor, conn, "3DT0", properties};
  my_io.BuildHierarchicalMatrix(32, 0, 1.e-3);
  std::cout << my_io.GetCompressionRatio()<<"\n";
  auto colpts = my_io.GetCollocationPoints();
  std::cout << colpts[4] <<"\n";
  ASSERT_TRUE((abs(colpts[0]-0.3333) < 1e-4 )  && (abs(colpts[1]-0.3333) < 1e-4)  && (abs(colpts[3]-0.6666) < 1e-4 )  && (abs(colpts[4]-0.6666) < 1e-4) );
}
/* -------------------------------------------------------------------------- */

TEST(bigwham_io_gen_3d, 3DT0_2) {
    // 2 Triangle for the square [0,1]X[0,1]
    int n_elts = 2;
    std::vector<double> coor(3*4  , 0.);
    coor[3]=1.;
    coor[7]=1.; coor[8]=0.;
    coor[9]=1.;coor[10]=1;
    std::vector<int> conn(n_elts * 3, 0);
    int k = 0;
    for (int i = 0; i < n_elts; i++) {
        conn[k] = i;
        conn[k + 1] = i + 1;
        conn[k + 2] = i + 2;
        k = k + 3;
    }
    std::vector<double> properties{1., 0.};

    BigWhamIOGen my_io{coor, conn, "3DT0", properties};
    my_io.BuildHierarchicalMatrix(32, 0, 1.e-3);
    std::cout << my_io.GetCompressionRatio()<<"\n";
    auto colpts = my_io.GetCollocationPoints();
    std::cout << colpts[4] <<"\n";
    std::vector<double> obspts(3  , 0.);
    il::Array<double> dd(3  , 1.);
    obspts[2]=10.;
    auto stress =my_io.ComputeStresses(obspts,dd.view());
    for (int i=0;i<6;i++){
        std::cout << "stress: " << stress[i] <<"\n";
    }

    ASSERT_TRUE(abs(stress[2]-0.000148955)<1.e-4 );
}
/* -------------------------------------------------------------------------- */
