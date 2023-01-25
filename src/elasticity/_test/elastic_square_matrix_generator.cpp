//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 24.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include <src/core/BoundaryElement.h>
#include <src/core/BEMesh.h>
#include <src/core/SquareMatrixGenerator.h>
#include <src/hmat/cluster/cluster.h>
#include <src/elasticity/2d/BIE_elastostatic_segment_0_impls.h>
#include <src/elasticity/2d/ElasticS3DP0_element.h>
#include <src/hmat/hmatrix/Hmat.h>

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
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,test);
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
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,test);
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
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,test);
    M.set_permutation(permutation);
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

//
TEST(SquareMatGen,segment_0_Hmat_1){
    ///  simple mesh
    int n_elts=88;
    il::Array2D<double> coor{n_elts+1,2,0.};

    il::Array2D<il::int_t> conn{n_elts,2};
    double h=0.123;
    for (int i=0;i<n_elts+1;i++) {
        coor(i,0)=i*h;
    }
    for (int i=0;i<n_elts;i++){
        conn(i,0)= i;
        conn(i,1)=i+1;
    }
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(coor,conn,seg0);
    il::Array2D<double> xcol=my_mesh.getCollocationPoints();
    bie::ElasticProperties elas(1,0.3);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  ker(elas,coor.size(1));
    il::Array<double> prop{1,1000.};
    ker.setKernelProperties(prop);

    il::int_t max_leaf_size=32;
    bie::Cluster cluster = bie::cluster(max_leaf_size, il::io, xcol);
    il::Array<il::int_t> permutation =cluster.permutation;
//
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,ker);
    M.set_permutation(permutation);
    bie::Hmat<double>  h_;
    h_.toHmat(M,cluster,xcol,1,1.e-3);

    ASSERT_TRUE( h_.isBuilt() );//h_.isBuilt()
}


TEST(SquareMatGen,segment_0_Hmat_2){
    ///  simple mesh
    int n_elts=1200;
    il::Array2D<double> coor{n_elts+1,2,0.};

    il::Array2D<il::int_t> conn{n_elts,2};
    double L=1.;double h=2.*L/n_elts;
    for (int i=0;i<n_elts+1;i++) {
        coor(i,0)=i*h-L;
    }
    for (int i=0;i<n_elts;i++){
        conn(i,0)= i;
        conn(i,1)=i+1;
    }
    bie::Segment<0> seg0;
    bie::BEMesh<bie::Segment<0>> my_mesh(coor,conn,seg0);
    il::Array2D<double> xcol=my_mesh.getCollocationPoints();
    bie::ElasticProperties elas(1,0.0);
    bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>  ker(elas,coor.size(1));
    il::Array<double> prop{1,1000.};
    ker.setKernelProperties(prop);

    il::int_t max_leaf_size=32;
    bie::Cluster cluster = bie::cluster(max_leaf_size, il::io, xcol);
    il::Array<il::int_t> permutation =cluster.permutation;
//
    bie::SquareMatrixGenerator<double,bie::Segment<0>,bie::BIE_elastostatic<bie::Segment<0>,bie::Segment<0>,bie::ElasticKernelType::H>> M(my_mesh,ker);
    M.set_permutation(permutation);
    bie::Hmat<double>  h_;
    h_.toHmat(M,cluster,xcol,3,1.e-3);
    //simple opening mode...
    il::Array<double> x{M.size(1),0.0},y{M.size(1),0.0};
    for(il::int_t i=0;i<M.sizeAsBlocks(0);i++){
        x[2*permutation[i]+1]=4.0*sqrt(L*L-xcol(i,0)*xcol(i,0) );
    }
    y=h_.matvec(x);
    il::Array<double> rel_err{M.sizeAsBlocks(0),0.};
    for (il::int_t i=0;i<M.sizeAsBlocks(0);i++){
        rel_err[i]=sqrt((y[2*i+1]-1.)*(y[2*i+1]-1.));
       //std::cout << "rel x: " << rel_err[i] << "\n";
    }
   // std::cout << "Linf rel error " << il::norm(rel_err,il::Norm::Linf) <<"\n";
    //std::cout << "L2 rel error " << il::norm(rel_err,il::Norm::L2) <<"\n";
    std::cout << "Mean rel error " << il::mean(rel_err) <<"\n";
    ASSERT_TRUE( il::mean(rel_err)<0.05 );//h_.isBuilt()
}