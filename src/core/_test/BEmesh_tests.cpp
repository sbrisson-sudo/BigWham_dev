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

#include <src/core/elements/Segment.h>

#include <src/core/BEMesh.h>

TEST(bemesh_seg,seg_0_1){
    int n_elts=4;
    il::Array2D<double> coor{n_elts+1,2,0.};
    il::Array2D<double> eltC{n_elts,2,0.};
    il::Array2D<il::int_t> conn{n_elts,2};
    double h=1;
    for (int i=0;i<n_elts+1;i++) {
        coor(i,0)=i*h;
    }
    for (int i=0;i<n_elts;i++){
            conn(i,0)= i;
            conn(i,1)=i+1;
            eltC(i,0)=(coor(i,0)+coor(i+1,0))/2.;
    }
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(coor,conn,seg0);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bool t= true;
    for (int i=0;i<n_elts;i++){
        t = t && (abs(test(i,0)-eltC(i,0))<1.e-6);
    }
    ASSERT_TRUE( t );
}


TEST(bemesh_seg,seg_0_2){
    int n_elts=6;
    il::Array2D<double> coor{n_elts+1,2,0.};
    il::Array2D<double> eltC{n_elts,2,0.};
    il::Array2D<il::int_t> conn{n_elts,2};
    double h=0.123;
    for (int i=0;i<n_elts+1;i++) {
        coor(i,1)=i*h;
        coor(i,0)=i*h;
    }
    for (int i=0;i<n_elts;i++){
        conn(i,0)= i;
        conn(i,1)=i+1;
        eltC(i,0)=(coor(i,0)+coor(i+1,0))/2.;
        eltC(i,1)=(coor(i,1)+coor(i+1,1))/2.;
    }
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(coor,conn,seg0);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bool t= true;
    for (int i=0;i<n_elts;i++){
        t = t && (abs(test(i,1)-eltC(i,1))<1.e-6) && (abs(test(i,0)-eltC(i,0))<1.e-6);
    }
    ASSERT_TRUE( t );
}
