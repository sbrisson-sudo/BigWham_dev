//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 25.01.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cstdlib>
#include <memory>
//
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <cmath>

#include "core/elastic_properties.h"

#include "elasticity/fullspace_iso_3d_triangle/elastic_3dT0_element.h"

#include "elasticity/fullspace_iso_3d_triangle/elastic_3dT0_element_af.h"

TEST(T0_ker,AngDislocation_1){
    // test angular dislocation displacement in angular dislocation system
    il::StaticArray2D<double,3,3> test_ad=bie::AngDisDisp(2.4,-3., 0.2, il::pi/2., 0.15);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_ad(i,0) << " | " << test_ad(i,1) << " | "<< test_ad(i,2) << "\n";
    }
    ASSERT_NEAR(test_ad(0,2),0.0528023,0.001);
}


TEST(T0_ker,TDSetupD_1){
    // test TDsetup
    il::StaticArray<double, 3> x_obs{il::value,{2.4,-3.0,0.2}};

    il::Array<double> TriVertex{3,0.};
    il::Array<double> SideVec{3,0.};
    SideVec[1]=1;SideVec[2]=1;

    il::StaticArray2D<double,3,3>  test_ad=bie::TDSetupD(x_obs,il::pi/2.+il::pi, 0.25,TriVertex,SideVec);

    for (il::int_t i=0;i<3;i++){
        std::cout   << test_ad(i,0) << " | " << test_ad(i,1) << " | "<< test_ad(i,2) << "\n";
    }
    ASSERT_NEAR(test_ad(1,1),0.0720,0.0001);
}


TEST(T0_ker,TDSetupD_2){
    // test tdSetupD 2
    il::StaticArray<double, 3> x_obs{il::value,{0.2,1./3,1./3}};

    il::Array<double> TriVertex{3,0.};
    il::Array<double> SideVec{3,0.};
    SideVec[1]=1;SideVec[2]=1;TriVertex[1]=1;

    il::StaticArray2D<double,3,3>  test_ad=bie::TDSetupD(x_obs,il::pi/2.+il::pi, 0.25,TriVertex,SideVec);

//    il::StaticArray2D<double,3,3> test_ad=bie::AngDisDisp(2.,-3., 0.2, il::pi/2., 0.15);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_ad(i,0) << " | " << test_ad(i,1) << " | "<< test_ad(i,2) << "\n";
    }
    ASSERT_NEAR(test_ad(1,1),0.0284,0.0001);
}

TEST(T0_ker,TDSetupD_3){
    // test tdSetupD 2
    il::StaticArray<double, 3> x_obs{il::value,{0.2,1./3,1./3}};

    il::Array<double> TriVertex{3,0.};
    il::Array<double> SideVec{3,0.};
    SideVec[1]=1;SideVec[2]=1;TriVertex[1]=0;

    il::StaticArray2D<double,3,3>  test_ad=bie::TDSetupD(x_obs,il::pi/2.+il::pi, 0.25,TriVertex,SideVec);

//    il::StaticArray2D<double,3,3> test_ad=bie::AngDisDisp(2.,-3., 0.2, il::pi/2., 0.15);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_ad(i,0) << " | " << test_ad(i,1) << " | "<< test_ad(i,2) << "\n";
    }
    ASSERT_NEAR(test_ad(1,1),0.015887,0.0001);
}


TEST(T0_ker,TDSetupD_4){
    // test tdSetupD 2
    il::StaticArray<double, 3> x_obs{il::value,{ 0.3,2.3, -0.4000}};

    il::Array<double> TriVertex{3,0.};
    il::Array<double> SideVec{3,0.};
    SideVec[1]=-1.0; // 0,-1,0
    TriVertex[2]=1;// 0,0,1
    std::cout << " p        " << TriVertex[0] << ","  << TriVertex[1] << "," << TriVertex[2] << "\n";
    std::cout << " SideVec  " << SideVec[0] << ","  << SideVec[1] << "," << SideVec[2] << "\n";
    std::cout << " x_obs    " << x_obs[0] << ","  << x_obs[1] << "," << x_obs[2] << "\n";
    il::StaticArray2D<double,3,3>  test_ad=bie::TDSetupD(x_obs,1.5708, 0.25,TriVertex,SideVec);
    std::cout << " ---\n";
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_ad(i,0) << " | " << test_ad(i,1) << " | "<< test_ad(i,2) << "\n";
    }
    ASSERT_NEAR(test_ad(1,1),-0.01032,0.0001);
}

//
//TEST(T0,disp_tdcs_1){
//    // displacement in TDCS at center of elememt
//    il::StaticArray<double, 3> x_obs{il::value,{0.,1./3.,0.5+1./3.}};
//    il::StaticArray2D<double,3,3> P{0.};
//    // make a TD  x=0,   P2 origin      P2-P3 as e_y
//    P(0,2)=1.0; // P1
//    P(2,1)=1.0;  P(2,2)=1.0; // P3
////    for (il::int_t i=0;i<3;i++){
////        x_obs[i]+=sqrt(std::numeric_limits<double>::epsilon());
////    }
//    il::StaticArray2D<double,3,3>  test_tdcs=bie::Displacements_TDCS(x_obs,P,0.25);
//
//    for (il::int_t i=0;i<3;i++){
//        std::cout   << test_tdcs(i,0) << " | " << test_tdcs(i,1) << " | "<< test_tdcs(i,2) << "\n";
//    }
//    ASSERT_NEAR(test_tdcs(1,1),0.5,0.0000000001);
//}


TEST(T0_ker,disp_efcs_1){
    // displacement in EFCS at center of elememt
    il::StaticArray<double, 3> x_obs{il::value,{1./3.,1./3.,0.}};
    il::StaticArray2D<double,3,3> P{0.};
    // unit triangle
    P(1,0)=1.0; // P2
    P(2,1)=1.0;  // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.25);

    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(1,1),0.5,0.0000000001);
}


TEST(T0_ker,disp_efcs_2){
    // displacement in EFCS  a given point for the unit triangle
    il::StaticArray<double, 3> x_obs{il::value,{1.4,2.3,0.3}};
    il::StaticArray2D<double,3,3> P{0.};
    P(1,0)=1.0; // P2
    P(2,1)=1.0;  // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.25);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(1,1),0.00189317,0.00000001);
}


TEST(T0_ker,disp_efcs_3){
    // displacement in EFCS  a given point for the unit triangle
    il::StaticArray<double, 3> x_obs{il::value,{-0.4,1.02,-0.12}};
    il::StaticArray2D<double,3,3> P{0.};
    P(1,0)=1.0; // P2
    P(2,1)=1.0;  // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.5);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(0,0),-0.0142032,0.000001);
}




// test a non unit triangle with normal e_3
TEST(T0_ker,disp_efcs_4){
    // displacement in EFCS  a given point for the unit triangle
    il::StaticArray<double, 3> x_obs{il::value,{0.4,1.02,-0.12}};
    il::StaticArray2D<double,3,3> P{0.};
    P(1,0)=1.2; // P2
    P(2,1)=1.3;  P(2,0)=0.5; // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.3);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(1,1),-0.181757,0.000001);
}


// test a non unit  triangle but still with normal = e3
TEST(T0_ker,disp_efcs_5){
    // displacement in EFCS at a given point for a triangle
    il::StaticArray<double, 3> x_obs{il::value,{0.4,1.02,-0.12}};
    il::StaticArray2D<double,3,3> P{0.};
    P(1,0)=1.2; P(1,1)=0.2;// P2
    P(2,1)=1.3;  P(2,0)=0.5;   // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.3);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(0,0),-0.146457,0.000001);
}

// test a non unit inclined triangle -> does not seem to work properly !
TEST(T0_ker,disp_efcs_6){
    // displacement in EFCS at a given point for an inclined triangle
    il::StaticArray<double, 3> x_obs{il::value,{0.4,1.02,-0.12}};
    il::StaticArray2D<double,3,3> P{0.};
    P(1,0)=1.2; // P2
    P(2,1)=1.3;  P(2,0)=0.5;  P(2,2)=0.5; // P3

    il::StaticArray2D<double,3,3>  test_efcs=bie::Displacements_EFCS(x_obs,P,0.3);
    for (il::int_t i=0;i<3;i++){
        std::cout   << test_efcs(i,0) << " | " << test_efcs(i,1) << " | "<< test_efcs(i,2) << "\n";
    }
    ASSERT_NEAR(test_efcs(0,0),-0.0507796,0.000001);
}


TEST(T0_ker,disp_af_1){
// unit triangle
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=1./3.;
    obsP[1]=1./3.;
    auto res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);

    for (il::int_t i=0;i<3;i++){
        std::cout   << res(i,0)  << " | " << res(i,1)  << " | "<< res(i,2)  << "\n";
    }
    ASSERT_TRUE(res(0,0) ==-0.5);
}


TEST(T0_ker,disp_af_2){

// unit triangle
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.2;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=3.;
    obsP[1]=-0.3;
    auto res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);
    bool test = true;

    for (il::int_t i=0;i<3;i++){
        std::cout   << res(i,0)  << " | " << res(i,1)  << " | "<< res(i,2)  << "\n";
    }
    ASSERT_TRUE(res(1,1)==0);
}


TEST(T0_ker,disp_af_3){

// unit triangle  with normal along e2
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.2;
    elVertex(2,0)=1.3;  elVertex(2,1)=0.0; elVertex(2,2)=1.0;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=3.; obsP[2]=2.;

    il::StaticArray2D<double,3,3> res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);

    for (il::int_t i=0;i<3;i++){
        std::cout   << res(i,0)  << " | " << res(i,1)  << " | "<< res(i,2)  << "\n";
    }
    ASSERT_TRUE(res(1,1)==0);
}


TEST(T0_ker,disp_af_4){

// inclined  triangle - all checks not implemented
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.2;
    elVertex(2,0)=1.3;  elVertex(2,1)=0.3; elVertex(2,2)=1.0;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=3.; obsP[2]=2.;

    il::StaticArray2D<double,3,3> res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);

    for (il::int_t i=0;i<3;i++){
        std::cout   << res(i,0)  << " | " << res(i,1)  << " | "<< res(i,2)  << "\n";
    }
    ASSERT_TRUE(true);
}


TEST(T0_ker,disp_1){
    // unit triangle - comparison with mma
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=4.;
    obsP[2]=0.4;
    for (il::int_t i=0;i<3;i++){
        obsP[i]+=sqrt(std::numeric_limits<double>::epsilon());
    }
    auto res = bie::DisplacementKernelT0(obsP,elVertex,0.25);

    bool test = true;

    il::StaticArray2D<double, 3, 3> mma_res;
    mma_res(0,0)=-0.000733982;mma_res(0,1)=0.0000512986; mma_res(0,2)=0.000896384;
    mma_res(1,0)=0.0000512986;mma_res(1,1)= -0.000113193;mma_res(1,2)= -0.0000759645;
    mma_res(2,0)=-0.00103584;mma_res(2,1)=0.0000870897;mma_res(2,2)= -0.000114662;
    //std::cout << " T0 displacement test 1\n";
    for (il::int_t i;i<3;i++){
        for (il::int_t j=0;j<3;j++){
            test = test && (std::abs(res(i,j)-mma_res(i,j)) < 1.e-5 );
            std::cout << " c++ implt:" << res(i,j) << " mma :" << mma_res(i,j) <<"\n";
        }
    }
    ASSERT_TRUE(test);
}


TEST(T0_ker,disp_1_1){

    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=4.;
    obsP[2]=0.4;

    auto res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);

    bool test = true;

    il::StaticArray2D<double, 3, 3> mma_res;
    mma_res(0,0)=-0.000733982;mma_res(0,1)=0.0000512986; mma_res(0,2)=0.000896384;
    mma_res(1,0)=0.0000512986;mma_res(1,1)= -0.000113193;mma_res(1,2)= -0.0000759645;
    mma_res(2,0)=-0.00103584;mma_res(2,1)=0.0000870897;mma_res(2,2)= -0.000114662;
    std::cout << " T0 displacement test 1\n";
    for (il::int_t i;i<3;i++){
        for (il::int_t j=0;j<3;j++){
            test = test && (std::abs(res(i,j)-mma_res(i,j)) < 1.e-5 );
            std::cout << " c++ implt:" << res(i,j) << " mma :" << mma_res(i,j) <<"\n";
        }
    }
    ASSERT_TRUE(test);
}

// THIS TEST DOES NOT PASS !!
// test highligthing the problematic singularities of Fata 2013 integration
TEST(T0_ker,disp_2){
    // test highligthing the problematic singularities of Fata 2013 integration
    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=45.;
    obsP[2]=0.;
    for (il::int_t i=0;i<3;i++){
        obsP[i]+=sqrt(std::numeric_limits<double>::epsilon());
    }
    auto res = bie::DisplacementKernelT0(obsP,elVertex,0.25);
    bool test = true;

    il::StaticArray2D<double, 3, 3> mma_res;
    mma_res(0,0)=-1.55237e-14; mma_res(0,1)= 9.85458e-17; mma_res(0,2)=6.64744e-6;
    mma_res(1,0)=9.85458e-17;mma_res(1,1)=-2.21892e-15;mma_res(1,2)= 4.93277e-8;
    mma_res(2,0)=-6.64744e-6;mma_res(2,1)= 4.93277e-8;mma_res(2,2)= -2.21783e-15;
    std::cout << " T0 displacement test 2\n";

    for (il::int_t i;i<3;i++){
        for (il::int_t j=0;j<3;j++){
            test = test && (std::abs(res(i,j)-mma_res(i,j)) < 1.e-5 );
            std::cout << " c++ implt:" << res(i,j) << " mma :" << mma_res(i,j) <<"\n";
        }
    }
    ASSERT_TRUE(test);
}


TEST(T0_ker,disp_2_1){

    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(1,0)=1.;
    elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{0.};
    obsP[0]=45.;
    obsP[2]=0.;

    auto res = bie::DisplacementKernelT0_af(obsP,elVertex,0.25);
    bool test = true;

    il::StaticArray2D<double, 3, 3> mma_res;
    mma_res(0,0)=-1.55237e-14; mma_res(0,1)= 9.85458e-17; mma_res(0,2)=6.64744e-6;
    mma_res(1,0)=9.85458e-17;mma_res(1,1)=-2.21892e-15;mma_res(1,2)= 4.93277e-8;
    mma_res(2,0)=-6.64744e-6;mma_res(2,1)= 4.93277e-8;mma_res(2,2)= -2.21783e-15;
    std::cout << " T0 displacement test 2\n";

    for (il::int_t i;i<3;i++){
        for (il::int_t j=0;j<3;j++){
            test = test && (std::abs(res(i,j)-mma_res(i,j)) < 1.e-5 );
            std::cout << " c++ implt:" << res(i,j) << " mma :" << mma_res(i,j) <<"\n";
        }
    }
    ASSERT_TRUE(test);
}


TEST(T0_ker,disp_3_1){

    il::StaticArray2D<double, 3, 3> elVertex{0.};
    elVertex(0,0)=-0.1;  elVertex(0,1)=-0.2;
    elVertex(1,0)=1.3; elVertex(1,1)=0.3;
    elVertex(2,0)=-0.4; elVertex(2,1)=1.;

    il::StaticArray<double, 3> obsP{il::value,{1.24, 0.02, -0.12}};

    auto res = bie::DisplacementKernelT0_af(obsP,elVertex,0.35);
    bool test = true;

    il::StaticArray2D<double, 3, 3> mma_res;
    mma_res(0,0)=0.0225481; mma_res(0,1)= -0.0160498; mma_res(0,2)=0.00954313;
    mma_res(1,0)=-0.0160498;mma_res(1,1)=0.0395082;mma_res(1,2)=-0.00663014;
    mma_res(2,0)=-0.0193029;mma_res(2,1)= 0.0301412;mma_res(2,2)=0.00993545;
    std::cout << " T0 displacement test 3\n";

    for (il::int_t i=0;i<3;i++){
        for (il::int_t j=0;j<3;j++){
            test = test && (std::abs(res(i,j)-mma_res(i,j)) < 1.e-5 );
            std::cout << i << " - "<< j << " c++ implt:" << res(i,j) << " mma :" << mma_res(i,j) <<"\n";
        }
    }
    ASSERT_TRUE(test);
}





