//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 13.12.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//

#include <gtest/gtest.h>

#include "core/oldies/DomainMesh.h"


TEST(domain_Mesh_Test,t1) {
 // test centroid.

  // background wellMesh 2 Quads
  il::Array2D<double> nodesB{6,2,0.};
  nodesB(0,0)=-100.;
  nodesB(0,1)=-100.;
  nodesB(1,0)=-100.;
  nodesB(1,1)=100.;
  nodesB(2,0)=0.;
  nodesB(2,1)=100.;
  nodesB(3,0)=0.;
  nodesB(3,1)=-100.;
  nodesB(4,0)=100.;
  nodesB(4,1)=-100.;
  nodesB(5,0)=100.;
  nodesB(5,1)=100.;
  il::Array2D<il::int_t> connB{2,4,0};
  connB(0,0)=0;
  connB(0,1)=1;
  connB(0,2)=2;
  connB(0,3)=3;
  connB(1,0)=2;
  connB(1,1)=3;
  connB(1,2)=4;
  connB(1,3)=5;

  il::Array<il::int_t> matidB{2,0};
  matidB[1]=1;

  bie::DomainMesh bckmesh(nodesB,connB,matidB);

  auto  kk1 = bckmesh.elementCentroid(0);
  auto  kk2 = bckmesh.elementCentroid(1);

  std::cout << kk1[1] << "\n";

  ASSERT_TRUE((kk1[0]==-50.)&&(kk1[1]==0.)&&(kk2[0]==50.)&&(kk2[1]==0.));
}

TEST(domain_Mesh_Test,t2) {

  // background wellMesh 2 Quads
  il::Array2D<double> nodesB{6,2,0.};
  nodesB(0,0)=-100.;
  nodesB(0,1)=-100.;
  nodesB(1,0)=-100.;
  nodesB(1,1)=100.;
  nodesB(2,0)=0.;
  nodesB(2,1)=100.;
  nodesB(3,0)=0.;
  nodesB(3,1)=-100.;
  nodesB(4,0)=100.;
  nodesB(4,1)=-100.;
  nodesB(5,0)=100.;
  nodesB(5,1)=100.;
  il::Array2D<il::int_t> connB{2,4,0};
  connB(0,0)=0;
  connB(0,1)=1;
  connB(0,2)=2;
  connB(0,3)=3;
  connB(1,0)=2;
  connB(1,1)=3;
  connB(1,2)=4;
  connB(1,3)=5;

  il::Array<il::int_t> matidB{2,0};
  matidB[1]=1;

  bie::DomainMesh bckmesh(nodesB,connB,matidB);

  il::StaticArray<double,2> myxt{0.};
  myxt[0] = 1.; myxt[1] = 0.;

  il::int_t kk = bckmesh.locate(myxt);

  ASSERT_EQ(kk, 1); // point (1,0) is in the 2nd element, i.e. 2-1 =1 in c

}


TEST(domain_Mesh_Test,t3) {

  // background wellMesh 2 Quads
  double far_away = 1.e12;
  il::Array2D<double> nodesB{6,2,0.};
  nodesB(0,0)=-far_away;
  nodesB(0,1)=-far_away;
  nodesB(1,0)=-far_away;
  nodesB(1,1)=far_away;
  nodesB(2,0)=0.;
  nodesB(2,1)=far_away;
  nodesB(3,0)=0.;
  nodesB(3,1)=-far_away;
  nodesB(4,0)=far_away;
  nodesB(4,1)=-far_away;
  nodesB(5,0)=far_away;
  nodesB(5,1)=far_away;
  il::Array2D<il::int_t> connB{2,4,0};
  connB(0,0)=0;
  connB(0,1)=1;
  connB(0,2)=2;
  connB(0,3)=3;
  connB(1,0)=3;
  connB(1,1)=2;
  connB(1,2)=5;
  connB(1,3)=4;

  il::Array<il::int_t> matidB{2,0};
  matidB[1]=1;

  bie::DomainMesh bckmesh(nodesB,connB,matidB);

  il::StaticArray<double,2> myxt{0.};
  myxt[0] = 1.e-12; myxt[1] = 0.1;

  il::int_t kk = bckmesh.locate(myxt);

  ASSERT_EQ(kk, 1); // point (1,0) is in the 2nd element, i.e. 2-1 =1 in c

}

TEST(domain_Mesh_Test,t4) {
  // test locate and matid
  // background wellMesh 2 Quads
  il::Array2D<double> nodesB{6,2,0.};
  nodesB(0,0)=-100.;
  nodesB(0,1)=-100.;
  nodesB(1,0)=-100.;
  nodesB(1,1)=100.;
  nodesB(2,0)=0.;
  nodesB(2,1)=100.;
  nodesB(3,0)=0.;
  nodesB(3,1)=-100.;
  nodesB(4,0)=100.;
  nodesB(4,1)=-100.;
  nodesB(5,0)=100.;
  nodesB(5,1)=100.;
  il::Array2D<il::int_t> connB{2,4,0};
  connB(0,0)=0;
  connB(0,1)=1;
  connB(0,2)=2;
  connB(0,3)=3;
  connB(1,0)=2;
  connB(1,1)=3;
  connB(1,2)=4;
  connB(1,3)=5;

  il::Array<il::int_t> matidB{2,0};
  matidB[1]=25;

  bie::DomainMesh bckmesh(nodesB,connB,matidB);

  il::StaticArray<double,2> myxt{};
  myxt[0] = 1.0; myxt[1] = 0.1;

  il::int_t kk = bckmesh.locate(myxt);

  ASSERT_EQ(matidB[kk], 25);

}
