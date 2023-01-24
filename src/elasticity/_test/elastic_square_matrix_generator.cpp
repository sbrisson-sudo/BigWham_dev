//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 24.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include <src/core/BoundaryElement.h>
#include <src/core/BEMesh.h>
#include <src/core/SquareMatrixGenerator.h>
#include <src/elasticity/2d/BIE_elastostatic_segment_0_impls.h>
#include <src/elasticity/2d/ElasticS3DP0_element.h>

TEST(SquareMatGen,segment_0_1){
     /// single element
    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    il::Array2D<il::int_t> conn{1,2,0};
    conn(0,1)=1;
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(xy,conn,seg0);
    bie::Segment<0> source;
//    source.setSegment(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<il::int_t> permutation{1,0};
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,permutation,test);
    ASSERT_TRUE(M.size(1)==2 && M.size(0)==2 );
}

TEST(SquareMatGen,segment_0_2){
    /// single element
    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    il::Array2D<il::int_t> conn{1,2,0};
    conn(0,1)=1;
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(xy,conn,seg0);
    bie::Segment<0> source;
//    source.setSegment(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<il::int_t> permutation{1,0};
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,permutation,test);
    ASSERT_TRUE(M.blockSize()==2 && M.sizeAsBlocks(0)==my_mesh.numberOfElts() );
}


TEST(SquareMatGen,segment_0_3){
    /// single element
    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    il::Array2D<il::int_t> conn{1,2,0};
    conn(0,1)=1;
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(xy,conn,seg0);
    bie::Segment<0> source;
//    source.setSegment(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<double> prop{1,1000.};
    test.setKernelProperties(prop);
    il::Array<il::int_t> permutation{1,0};
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,permutation,test);
    il::Array2D<double> A{M.size(0),M.size(1),0.0};
    il::Array2DEdit<double> v=A.Edit();

    M.set(0,0,il::io,v);
    // check with known values of entry for that case
// 0.34979115367667662
//0.34979125861394933
//    std::cout << v(0,0) << "-" << v(0,1) <<"\n";
 //   std::cout << v(1,0) << "-" << v(1,1) <<"\n";
    ASSERT_TRUE( (abs(v(0,0)-0.34979115367667662)<1.e-12) && (abs(v(0,1)-0.0)<1.e-12) && (abs(v(1,0)-0.0)<1.e-12) && (abs(v(1,1)-0.34979125861394933)<1.e-12) );
}
