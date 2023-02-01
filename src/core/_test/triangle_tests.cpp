//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#pragma once
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include <src/core/elements/Triangle.h>

TEST(triangle,triangle_0_1){

    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);

    ASSERT_TRUE(tri0.getSpatialDimension()==3);
}

TEST(triangle,triangle_0_2){
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    ASSERT_TRUE(tri0.getNumberOfCollocationPoints()==1);
}

TEST(triangle,triangle_0_3){
    // check rotation matrix
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto R = tri0.rotationMatrix();
    ASSERT_TRUE(R(0,0)==1 && R(1,1)==1 && R(2,2)==1 && R(0,1)==0 && R(0,2)==0& R(1,0)==0 && R(1,2)==0 && R(2,0)==0 && R(2,1)==0);
}

TEST(triangle,triangle_0_4){
    // check centroid
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getCentroid();
    ASSERT_TRUE(x[0]==1./3&& x[1]==1./3 && x[2]==0.);
}


TEST(triangle,triangle_0_5){
    // check nodes
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getNodes();
    ASSERT_TRUE(abs(x(0,0)-1./3)<1.e-12 && abs(x(0,1)-1./3)<1.e-12 && abs(x(0,2))<1.e-12);
}


TEST(triangle,triangle_0_6){
    // check vertices
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getVertices();
    ASSERT_TRUE(x(0,0)==xyz(0,0) && x(0,1)==xyz(0,1) && x(0,2)==xyz(0,2) );
}



TEST(triangle,triangle_0_7){
    // check normal
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getNormal();
    ASSERT_TRUE(x[0]==0 && x[1]==0 && x[2]==1 );
}


TEST(triangle,triangle_0_8){
    // check tangent 1
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getTangent_1();
    ASSERT_TRUE(x[0]==1 && x[1]==0 && x[2]==0 );
}


TEST(triangle,triangle_0_9){
    // check tangent 2
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<0> tri0;
    tri0.setElement(xyz);
    auto x= tri0.getTangent_2();
    ASSERT_TRUE(x[0]==0 && x[1]==1 && x[2]==0 );
}






TEST(triangle,triangle_2_0){
    // check # nodes
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<2> tri2;
    tri2.setElement(xyz);
    auto x= tri2.getNodes();
    ASSERT_TRUE(tri2.getNumberOfNodes()==6);
}


TEST(triangle,triangle_2_1){
    //  check dimension
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<2> tri2;
    tri2.setElement(xyz);
    auto x= tri2.getNodes();
    std::cout  << x(0,0) ;
    ASSERT_TRUE(tri2.getSpatialDimension()==3);
}


TEST(triangle,triangle_2_2){
    // check # collocation points
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<2> tri2;
    tri2.setElement(xyz);
    auto x= tri2.getNodes();
    std::cout  << x(0,0) ;
    ASSERT_TRUE(tri2.getNumberOfCollocationPoints()==6);
}


TEST(triangle,triangle_2_3){
    // check rotation matrix
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<2> tri2;
    tri2.setElement(xyz);
    auto R = tri2.rotationMatrix();
    ASSERT_TRUE(R(0,0)==1 && R(1,1)==1 && R(2,2)==1 && R(0,1)==0 && R(0,2)==0& R(1,0)==0 && R(1,2)==0 && R(2,0)==0 && R(2,1)==0);
}


TEST(triangle,triangle_2_4){
    // check collocation points location
    il::Array2D<double> xyz{3,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    bie::Triangle<2> tri2;
    tri2.setElement(xyz);
    auto xcol = tri2.getCollocationPoints();
    ASSERT_TRUE(xcol.size(0)==6 && xcol.size(1)==3);
}

