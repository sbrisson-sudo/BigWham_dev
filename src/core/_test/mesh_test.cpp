//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 24.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include "core/Mesh2D.h"

/// TEST 1
TEST(Mesh, test_ribbon) {
  //  a very simple wellMesh with 4 elements  (0,1,2,3)
  // ensure the ribbon elements are 1 and 2

  // create the wellMesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  //  // create a basic 1D wellMesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  bie::Mesh mesh(xy, myconn, 0);
  //
  auto ribbon = mesh.getRibbonElements();
  //
  //  std::cout << "ribbon size: " << ribbon.size() << " "<<ribbon[0] << " " <<
  //  ribbon[1];
  //

  ASSERT_TRUE((ribbon[0] == 1) && (ribbon[1] == 2));  //&& (ribbon[1]==2)
}

/// TEST 2
TEST(Mesh, test_mesh_object) {
  //  a very simple Mesh with 4 elements  (0,1,2,3)
  // ensure the mesh class works correctly

  // create the Mesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::int_t e = 0;

  bie::Mesh mesh(xy, myconn, 0);

  ASSERT_TRUE((mesh.numberOfElts() == nelts) &&
              (mesh.numberOfElts() == myconn.size(0)) &&
              (mesh.numberOfNodes() == xy.size(0)) &&
              (mesh.interpolationOrder() == 0) && (mesh.tipElts(0) == 0) &&
              (mesh.tipElts(1) == 3) && (mesh.tipNodes(0) == 0) &&
              (mesh.tipNodes(1) == 4) && (mesh.eltSize(e) == 1.) &&
              (mesh.allEltSize()[0] == 1.) && (mesh.numberOfFractures() == 1));
}

// test add elements
TEST(Mesh, add_elt) {
  // test the addition of 2 elements at 90 degree from the left tip

  // create the Mesh.
  il::int_t nelts = 2;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::int_t e = 0;

  bie::Mesh mesh(xy, myconn, 0);

  mesh.addNTipElts(1, 2, 2, il::pi / 2.);

  ASSERT_TRUE(
      (mesh.numberOfElts() == nelts + 2) &&
      (mesh.coordinates(nelts + 1, 1) == 1) &&
      (mesh.coordinates(nelts + 2, 1) == 2) &&
      (mesh.coordinates(nelts + 2, 0) == mesh.coordinates(nelts + 1, 0)) &&
      (mesh.fracid(nelts + 1) == mesh.fracid(nelts)) &&
      (mesh.fracid(nelts + 1) == 0));
}

// test add elements
TEST(Mesh, add_elt_both_tip) {
  // test the addition of 2 elements at 90 degree from the left tip
   // check the addition of elts in both tip and check
  // consistency of tip position independent of the order of addition
  // create the Mesh.
  il::int_t nelts = 2;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::int_t e = 0;

  bie::Mesh mesh(xy, myconn, 0);
  bie::Mesh mesh2(xy, myconn, 0);

  mesh.addNTipElts(1, 2, 2, il::pi / 2.);
  mesh.addNTipElts(0, 0, 2, il::pi / 2.);

  mesh2.addNTipElts(0, 0, 2, il::pi / 2.);
  mesh2.addNTipElts(1, 2, 2, il::pi / 2.);
//  mesh2.tipElts(0) == mesh.tipElts(1) &&
//      mesh2.tipElts(1) == mesh.tipElts(0) &&
  ASSERT_TRUE(
          mesh2.coordinates(mesh2.tipElts(0),0)==mesh.coordinates(mesh.tipElts(0),0) &&
      mesh2.coordinates(mesh2.tipElts(0),1)==mesh.coordinates(mesh.tipElts(0),1) &&
              mesh2.coordinates(mesh2.tipElts(1),0)==mesh.coordinates(mesh.tipElts(1),0) &&
              mesh2.coordinates(mesh2.tipElts(1),1)==mesh.coordinates(mesh.tipElts(1),1)
  );
}



TEST(Mesh, elt_from_dof) {
  // test the method that get the element for which the  dof # given belongs
  il::int_t nelts = 2;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  bie::Mesh mesh(xy, myconn, 0);

  ASSERT_TRUE(mesh.eltFromdofDD(3)==1);

}

TEST(Mesh, colPts_p1) {
  // test the method getting the coordinates of all collocation points coor in the mesh
  // test for p1
  il::int_t nelts = 2;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::Array2D<double> col_pts{2*nelts, 2, 0.};
  il::int_t j=0;
  double xmean;
  for (il::int_t i = 0; i < nelts ; ++i) {
    xmean = 0.5*(xy(myconn(i,0),0)+xy(myconn(i,1),0));
    col_pts(j,0)=xmean-0.5/sqrt(2.);
    j++;
    col_pts(j,0)=xmean+0.5/sqrt(2.);
    j++;

  }

  bie::Mesh mesh(xy, myconn, 1);
  il::Array2D<double> col_pts_2 =mesh.getCollocationPoints();

  std::cout << col_pts(0,0) <<  " " << col_pts(0,1)  <<"\n";
  std::cout << col_pts_2(0,0) <<  " " << col_pts_2(0,1)   <<"\n";

  ASSERT_TRUE( (col_pts(0,0)==col_pts_2(0,0) ) && (col_pts(1,0)==col_pts_2(1,0))
                   && (col_pts(2,0)==col_pts_2(2,0) ) && (col_pts(3,0)==col_pts_2(3,0) )
  );

}


TEST(Mesh, col_dof_p1) {
  // test the method getting the dof map by  collocation points  in the mesh
  // test for p1
  il::int_t nelts = 2;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::Array2D<il::int_t > col_dof{2*nelts, 2, 0};
  il::int_t jk=0;
  for (il::int_t i = 0; i < nelts ; ++i) {
    // case P1 2 col per elements
        for (il::int_t j=0;j<2;j++){
          col_dof(i*2+j,0)=jk;
          col_dof(i*2+j,1)=jk+1;
          jk=jk+2;
        }
  }

  bie::Mesh mesh(xy, myconn, 1); // P1
  il::Array2D<il::int_t > col_dof_2 =mesh.getCollocationPointsDDDofs();

  for (il::int_t i=0;i<4;i++) {
    std::cout << col_dof(i, 0) << " " << col_dof(i, 1) << "\n";
    std::cout << col_dof_2(i, 0) << " " << col_dof_2(i, 1) << "\n";
  }

  ASSERT_TRUE( (col_dof(0,0)==col_dof_2(0,0) ) && (col_dof(1,0)==col_dof_2(1,0))
                   && (col_dof(2,0)==col_dof_2(2,0) ) && (col_dof(3,0)==col_dof_2(3,0) )  );

}



// todo Mesh tests for FractureID, MatID related as well as tip ordering
// related.
