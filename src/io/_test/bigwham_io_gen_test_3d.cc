//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 31.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "hmat/hmatrix/hmat.h"
#include "io/bigwham_io.h"

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
  int n_openMP_threads = 4;
  int n_GPUs = -1;
  bool verbose = true;
  bool homogeneous_size_pattern = false;
  bool useCuda = false;
  BigWhamIO my_io{coor, conn, "3DT0-H", properties, n_openMP_threads, n_GPUs, verbose, homogeneous_size_pattern, useCuda};
  my_io.BuildHierarchicalMatrix(32, 0, 1.e-3);
  std::cout << my_io.GetCompressionRatio()<<"\n";
  auto colpts = my_io.GetCollocationPoints();
  std::cout << colpts[4] <<"\n";
  ASSERT_TRUE((abs(colpts[0]-0.3333) < 1e-4 )  && (abs(colpts[1]-0.3333) < 1e-4)  && (abs(colpts[3]-0.6666) < 1e-4 )  && (abs(colpts[4]-0.6666) < 1e-4) );
}
/* -------------------------------------------------------------------------- */

TEST(bigwham_io_gen_3d, 3DT0_2) {
    // 2 Triangles for the square [0,1]X[0,1]
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
    int n_openMP_threads = 4;
    int n_GPUs = -1;
    bool verbose = true;
    bool homogeneous_size_pattern = false;
    bool useCuda = false;
    BigWhamIO my_io{coor, conn, "3DT0-H", properties, n_openMP_threads, n_GPUs, verbose, homogeneous_size_pattern, useCuda};
    my_io.BuildHierarchicalMatrix(32, 0, 1.e-3);
    std::cout << my_io.GetCompressionRatio()<<"\n";
    auto colpts = my_io.GetCollocationPoints();
    std::vector<double> obspts(3  , 0.);
    il::Array<double> dd(6  , 1.);
    obspts[2]=10.;
    auto stress =my_io.ComputeStresses(obspts,dd.view());
    for (int i=0;i<6;i++){
        std::cout << "stress: " << stress[i] <<"\n";
    }

    ASSERT_TRUE(abs(stress[2]-0.000299639)<1.e-4 );
}
/* -------------------------------------------------------------------------- */



TEST(bigwham_io_gen_3d, 3DT0_3) {
    // 2 Triangles for the square [0,1]X[0,1]
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
    int n_openMP_threads = 4;
    int n_GPUs = -1;
    bool verbose = true;
    bool homogeneous_size_pattern = false;
    bool useCuda = false;

    BigWhamIO my_io{coor, conn, "3DT0-H", properties, n_openMP_threads, n_GPUs, verbose, homogeneous_size_pattern, useCuda};
    my_io.BuildHierarchicalMatrix(32, 0, 1.e-3);
    std::cout << "Compression ratio:" << my_io.GetCompressionRatio()<<"\n";
    auto colpts = my_io.GetCollocationPoints();
    std::vector<double> obspts(3  , 0.);
    il::Array<double> dd(6  , 1.);
    obspts[2]=10.;obspts[0]=0.3;
    obspts[1]=-2.44;

    std::cout <<" dof dim " << my_io.dof_dimension() <<"\n";

    il::Array<double> displ =my_io.ComputeDisplacements(obspts,dd.view());
    for (int i=0;i<displ.size();i++){
        std::cout << "displacement: " << displ[i] <<"\n";
    }

    ASSERT_TRUE(abs(displ[0]+0.000158716)<1.e-4 );
}
/* -------------------------------------------------------------------------- */



TEST(bigwham_io_gen_3d, 3DR0_1) {
    // 1 Rectangle for the square [0,1]X[0,1]
    int n_elts =  1;
    int n_vertex = 4;
    std::vector<double> coor(3*4  , 0.);
    double h=0.5;
    coor[0]=-0.5; coor[1]=-0.5;
    coor[3]=-0.5;coor[4]=0.5;
    coor[6]=0.5;coor[7]=0.5;
    coor[9]=0.5; coor[10]=-0.5;
    std::vector<int> conn(n_elts * 4, 0);
    int k = 0;
    for (int i = 0; i < n_elts; i++) {
        conn[k] = i;
        conn[k + 1] = i + 1;
        conn[k + 2] = i + 2;
        conn[k + 3] = i + 3;
        k = k + 4;
    }
    std::vector<double> properties{1., 0.};
    int n_openMP_threads = 4;
    int n_GPUs = -1;
    bool verbose = true;
    bool homogeneous_size_pattern = false;
    bool useCuda = false;
    BigWhamIO my_io{coor, conn, "3DR0-H", properties, n_openMP_threads, n_GPUs, verbose, homogeneous_size_pattern, useCuda};
    std::cout << my_io.kernel_name() <<"\n";
    my_io.BuildHierarchicalMatrix(32, 0.0, 1.e-3);
    std::cout << my_io.GetCompressionRatio()<<"\n";
    auto colpts = my_io.GetCollocationPoints();
    std::cout << colpts[0] <<"\n";
    ASSERT_TRUE((abs(colpts[0]-0.) < 1e-4 )  && (abs(colpts[1]-0.) < 1e-4)  && (abs(colpts[3]-0.) < 1e-4 ) );
}
/* -------------------------------------------------------------------------- */

