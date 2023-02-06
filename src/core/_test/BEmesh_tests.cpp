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

#include "core/BEMesh.h"
#include "core/elements/Segment.h"
#include "core/elements/Triangle.h"
#include "core/elements/Rectangle.h"

//// 2D mesh test
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
    bie::BEMesh<bie::Segment<0>> my_mesh(coor, conn);
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
    bie::BEMesh<bie::Segment<0>> my_mesh(coor, conn);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bool t= true;
    for (int i=0;i<n_elts;i++){
        t = t && (abs(test(i,1)-eltC(i,1))<1.e-6) && (abs(test(i,0)-eltC(i,0))<1.e-6);
    }
    ASSERT_TRUE( t );
}


TEST(bemesh_seg,seg_1_1){
    int n_elts=4;
    il::Array2D<double> coor{n_elts+1,2,0.};
    il::Array2D<double> eltC{n_elts,2,0.};
    il::Array2D<il::int_t> conn{n_elts,2};
    double h=7.23;
    for (int i=0;i<n_elts+1;i++) {
        coor(i,0)=i*h;
    }
    for (int i=0;i<n_elts;i++){
        conn(i,0)= i;
        conn(i,1)=i+1;
        eltC(i,0)=(coor(i,0)+coor(i+1,0))/2.;
    }
    bie::BEMesh<bie::Segment<1>> my_mesh(coor, conn);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bool t= true;
    for (int i=0;i<n_elts;i++){ // left col points
        t = t && (abs(test(2*i,0)-eltC(i,0)+(h/2.)/sqrt(2))<1.e-6);
    }
    for (int i=0;i<n_elts;i++){ // right col points
        t = t && (abs(test(2*i+1,0)-eltC(i,0)-(h/2.)/sqrt(2))<1.e-6);
    }
    ASSERT_TRUE( t && (test.size(0)==n_elts*2) );//
}

////// triangular Mesh tests

TEST(bemesh_triangle,triangle_0_1){
    // mesh consisting of 2 triangles on the e_3 plane
    il::Array2D<double> xyz{4,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    xyz(3,1)=1.;xyz(3,0)=1.; xyz(3,2)=0.;
    il::Array2D<il::int_t> conn{2,3,0};
    conn(0,1)=1;conn(0,2)=2;
    conn(1,0)=2;conn(1,1)=1;conn(1,2)=3;
    bie::BEMesh<bie::Triangle<0>> my_mesh(xyz, conn);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    ASSERT_TRUE(test.size(0)==my_mesh.numberCollocationPoints()  &&  abs(test(0,2) ) <1.e-10 &&  abs(test(0,1) - 1./3.) <1.e-10 &&
                        abs(test(0,0) - 1./3.) <1.e-10 && abs(test(1,0) - 2./3.) <1.e-10 && abs(test(1,1) - 2./3.) <1.e-10 && abs(test(1,2) ) <1.e-10 );
}

TEST(bemesh_triangle,triangle_0_2){
    // mesh consisting of 2 triangles at 90 deg from one another
    il::Array2D<double> xyz{4,3,0.0};
    xyz(1,0)=1.;
    xyz(2,1)=1.;
    xyz(3,1)=1.;xyz(3,0)=0.; xyz(3,2)=1.; // at 90 degree
    il::Array2D<il::int_t> conn{2,3,0};
    conn(0,1)=1;conn(0,2)=2;
    conn(1,0)=2;conn(1,1)=1;conn(1,2)=3;
    bie::BEMesh<bie::Triangle<0>> my_mesh(xyz, conn);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bie::Triangle<0> tri_a;
    il::Array2D<double> xv=my_mesh.getVertices(0);
    tri_a.setElement(xv);
    xv=my_mesh.getVertices(1);
    bie::Triangle<0> tri_b;tri_b.setElement(xv);
    auto ndots= il::dot(tri_a.getNormal(),tri_b.getNormal());
    std::cout << " normal elt 1::" << tri_a.getNormal()[0] <<"-" << tri_a.getNormal()[1]<<"-"<<  tri_a.getNormal()[2] <<"\n";
    std::cout << " normal elt 2::" << tri_b.getNormal()[0] <<"-" <<tri_b.getNormal()[1] <<"-" <<tri_b.getNormal()[2] << "\n";
    auto ndots_1_1= il::dot(tri_a.getNormal(),tri_a.getTangent_1());
    auto ndots_1_2= il::dot(tri_a.getNormal(),tri_a.getTangent_2());
    auto ndots_2_1= il::dot(tri_b.getNormal(),tri_b.getTangent_1());
    auto ndots_2_2= il::dot(tri_b.getNormal(),tri_b.getTangent_2());
    ASSERT_TRUE(ndots==0 && ndots_1_1==0 && ndots_1_2==0 && ndots_2_1==0 && ndots_2_2==0);
}

// Rectangular mesh

TEST(bemesh_rectangle,rect_0_1){
    // cartesian mesh with _x*n_y rectangles on the same plane
    il::int_t n_x=3,n_y=1;il::int_t nelts=n_x*n_y;
    double hx=1.;double hy=2.;
    il::Array2D<double> coor{(n_x+1)*(n_y+1),3,0.};
    for (il::int_t k2=0;k2<n_y+1;k2++){
        for (il::int_t k=0;k<n_x+1;k++){
            coor(k+k2*(n_x+1),0)=hx*k;
            coor(k+k2*(n_x+1),1)=hy*k2;
        }
    }
    il::Array2D<il::int_t> conn{nelts,4,0};
    il::int_t c=0;
    for (il::int_t k2=0;k2<n_y;k2++){
        for (il::int_t k=0;k<n_x;k++){
            conn(c,0)=k+k2*n_y;conn(c,1)=k+1+k2*n_y;
            conn(c,2)=(k+1)+k2*n_y+n_x+1;conn(c,3)=(k)+k2*n_y+n_x+1;
            c++;
        }
    }
    for (il::int_t c=0;c<nelts;c++){
        std::cout << " elt #" << c<< " conn - " << conn(c,0) <<"-"<< conn(c,1) <<"-" <<conn(c,2) << "-" << conn(c,3) <<"\n";
    }
    bie::BEMesh<bie::Rectangle<0>> my_mesh(coor, conn);
    il::Array2D<double> test=my_mesh.getCollocationPoints();
    bie::Rectangle<0> rec_a;
    il::Array2D<double> xv=my_mesh.getVertices(0);

    rec_a.setElement(xv);
    auto aux=il::dot(rec_a.getTangent_1(),rec_a.getNormal());
    auto aux2=il::dot(rec_a.getTangent_1(),rec_a.getTangent_2());
    std::cout <<  rec_a.area() <<"\n";
    ASSERT_TRUE(my_mesh.numberOfElts()==nelts && aux==0 && aux2==0 && (rec_a.area()==2.0) && (test.size(0)==my_mesh.numberCollocationPoints()));
}


