//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include <src/core/BoundaryElement.h>
#include <src/elasticity/2d/BIE_elastostatic_segment_0_impls.h>

#include <src/elasticity/2d/ElasticS3DP0_element.h>

TEST(SP3D,test_seg_0_dof_dim){

    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    ASSERT_TRUE(test.getDofDimension() == 2);
}

TEST(SP3D,test_seg_0_dim){
    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    ASSERT_TRUE(test.getSpatialDimension()==2);
}



TEST(SP3D,test_seg_0_self){

    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<double> prop{1,1000.};
    test.setKernelProperties(prop);
    std::vector<double> test_self = test.influence(source,0,source,0);
    std::cout << "test self effect " << "\n";
    for (int i=0;i<4;i++){
        std::cout << test_self[i]  <<"\n";
    }
// old way // would have to be removed at some point !
    il::StaticArray2D<double,2,2> xys{0.};
    xys(1,0)=1.0;
    bie::SegmentData  auxi(xys,0);
    il::StaticArray2D<double,2,2> stnl = bie::normal_shear_stress_kernel_s3d_dp0_dd_nodal(auxi,auxi,0,0,elas,prop[0]);
    for (int i=0;i<2;i++){
        std::cout << stnl(i,0) <<"-"  <<stnl(i,1)  <<"\n";
    }
    ASSERT_TRUE( abs(stnl(0,0)-test_self[0])<1.e-5 && abs(stnl(1,0)-test_self[1])<1.e-5 &&  abs(stnl(0,1)-test_self[2])<1.e-5 &&  abs(stnl(1,1)-test_self[3])<1.e-5)  ;
}


TEST(SP3D,test_seg_0_1){
    il::Array2D<double> xy{2,2,0.};
    xy(0,1)=2.4; xy(1,0)=3.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    std::cout << "test on inclined elt " << "\n";
    il::Array2D<double> xy_r{2,2,0.};
    xy_r(0,0)=1.0; xy_r(1,0)=5.0; xy_r(0,1)=1.0;xy_r(1,1)=0.0;
    bie::Segment<0> receiver;
    receiver.setElement(xy_r);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<double> prop{1,1000.};
    test.setKernelProperties(prop);
    std::vector<double> test_self = test.influence(source,0,receiver,0);
// old way ..... // would have to be removed at some point !
    il::StaticArray2D<double,2,2> xys{0.},xys_r{0.};
    for (int i=0;i<2;i++){
        xys(i,0)=xy(i,0); xys(i,1)=xy(i,1);
        xys_r(i,0)=xy_r(i,0); xys_r(i,1)=xy_r(i,1);
    }
    bie::SegmentData  auxi_s(xys,0);
    bie::SegmentData  auxi_r(xys_r,0);
    il::StaticArray2D<double,2,2> stnl = bie::normal_shear_stress_kernel_s3d_dp0_dd_nodal(auxi_s,auxi_r,0,0,elas,prop[0]);
//    for (int i=0;i<2;i++){
//        std::cout << stnl(i,0) <<"-"  <<stnl(i,1)  <<"\n";
//    }
    ASSERT_TRUE( abs(stnl(0,0)-test_self[0])<1.e-5 && abs(stnl(1,0)-test_self[1])<1.e-5 &&  abs(stnl(0,1)-test_self[2])<1.e-5 &&  abs(stnl(1,1)-test_self[3])<1.e-5)  ;
}