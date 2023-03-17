//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "core/BoundaryElement.h"
#include "elasticity/FsIso2dSegment/BIE_elastostatic_segment_0_impls.h"
//#include "elasticity/FsIsoSp3dSegment/ElasticS3DP0_element.h"
#include <elasticity/FsIsoSp3dSegment/BieElastostaticSp3d.h>

TEST(SP3D,test_seg_0_dof_dim){

    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BieElastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    ASSERT_TRUE(test.getDofDimension() == 2);
}

TEST(SP3D,test_seg_0_dim){
    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BieElastostaticSp3d<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    ASSERT_TRUE(test.getSpatialDimension()==2);
}



TEST(SP3D,test_seg_0_self){

    il::Array2D<double> xy{2,2,0.};
    xy(1,0)=1.0;
    bie::Segment<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1,0.3);
    bie::BieElastostaticSp3d<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<double> prop{1,1000.};
    test.setKernelProperties(prop);
    std::vector<double> test_self = test.influence(source,0,source,0);
    std::cout << "test self effect " << "\n";
    for (int i=0;i<4;i++){
        std::cout << test_self[i]  << "\n";
    }
// old way // would have to be removed at some point !
    il::StaticArray2D<double,2,2> xys{0.};
    xys(1,0)=1.0;
    bie::SegmentData  auxi(xys,0);
    il::StaticArray2D<double,2,2> stnl = bie::normal_shear_stress_kernel_s3d_dp0_dd_nodal(auxi,auxi,0,0,elas,prop[0]);
    for (int i=0;i<2;i++){
        std::cout << stnl(i,0) <<"-"  <<stnl(i,1)  <<"\n";
    }
    double epsilon=1.e-6;
    ASSERT_TRUE( abs(stnl(0,0)-test_self[0])<epsilon && abs(stnl(1,0)-test_self[1])<epsilon &&  abs(stnl(0,1)-test_self[2])<epsilon &&  abs(stnl(1,1)-test_self[3])<epsilon)  ;
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
    bie::BieElastostaticSp3d<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  test(elas,xy.size(1));
    il::Array<double> prop{1,1000.};
    test.setKernelProperties(prop);
    std::vector<double> test_self = test.influence(source,0,receiver,0);
//    for (int i=0;i<3;i++){
//        std::cout << test_self[i] <<"-"  "\n";
//    }
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
    double eps=3.e-5;
    ASSERT_TRUE( abs(stnl(0,0)-test_self[0])<eps && abs(stnl(1,0)-test_self[1])<eps &&  abs(stnl(0,1)-test_self[2])<eps &&  abs(stnl(1,1)-test_self[3])<eps)  ;
}
