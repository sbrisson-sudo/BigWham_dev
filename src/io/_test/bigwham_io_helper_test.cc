//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.12.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "hmat/hmatrix/hmat.h"
#include "io/bigwham_io_helper.h"

#include "core/be_mesh.h"
#include "elements/boundary_element.h"
#include "elements/point.h"


TEST(mesh_from_vect,test_pts_2d_u){
    int n_elts = 4;
    std::vector<double> coor(n_elts*2);
    std::vector<int> conn(n_elts,0);
    std::cout << "creating mesh obj \n";
    double h = 1;
    int k=0;
    for (int i = 0; i < n_elts ; i++) {
        coor[k] = i * h;
        coor[k+1] = i * 2* h;
        conn[i]=i;
        k=k+2;
    }
    std::cout << "creating mesh obj \n";
    using EltType = bie::Point<2>;
    std::unique_ptr<bie::Mesh> mesh_obs=bie::CreateUniqueMeshFromVect<EltType>(2, 1, coor, conn);
    bool t = true;
    auto colpts =  mesh_obs->collocation_points();
    std::cout << "n elts :: " << mesh_obs->num_elements() <<"\n";
    for (int i = 0; i < n_elts; i++) {
         auto receiver_element = mesh_obs->GetElement(i);
         auto colpts_e = receiver_element->collocation_points();
         t= t && ( colpts(i,0) == i * h )  && (colpts(i,1)==i * 2 * h);
         t = t && ( colpts_e(0, 0) ==i * h ) && ( colpts_e(0, 1) ==i * 2* h)  ;
    }
    ASSERT_TRUE(t);

}


TEST(mesh_from_vect,test_pts_2d_s){
    int n_elts = 4;
    std::vector<double> coor(n_elts*2);
    std::vector<int> conn(n_elts,0);
    std::cout << "creating mesh obj \n";
    double h = 1;
    int k=0;
    for (int i = 0; i < n_elts ; i++) {
        coor[k] = i * h;
        coor[k+1] = i * 2* h;
        conn[i]=i;
        k=k+2;
    }
    std::cout << "creating mesh obj \n";
    using EltType = bie::Point<2>;
    std::shared_ptr<bie::Mesh> mesh_obs= bie::CreateMeshFromVect<EltType>(2, 1, coor, conn);
    bool t = true;
    auto colpts =  mesh_obs->collocation_points();
    std::cout << "n elts :: " << mesh_obs->num_elements() <<"\n";
    for (int i = 0; i < n_elts; i++) {
        auto receiver_element = mesh_obs->GetElement(i);
        auto colpts_e = receiver_element->collocation_points();
        t= t && ( colpts(i,0) == i * h )  && (colpts(i,1)==i * 2 * h);
        t = t && ( colpts_e(0, 0) ==i * h ) && ( colpts_e(0, 1) ==i * 2* h)  ;
    }
    ASSERT_TRUE(t);
}