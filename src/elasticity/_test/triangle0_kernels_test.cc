//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <cstdlib>
#include <memory>

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <cmath>

#include "core/elastic_properties.h"
#include "elasticity/bie_elastostatic.h"
#include "elasticity/fullspace_iso_3d_triangle/bie_elastostatic_triangle_0_influence.h"
#include "elements/triangle.h"
#include "elements/point.h"


TEST(Triangle0, test_H_1) {
  // test self-effect triangle 0
  il::Array2D<double> xy{3, 3, 0.};
  xy(1, 0) = 1;
  xy(2, 1) = 1.0;
  bie::Triangle<0> source;
  source.SetElement(xy);
  bie::ElasticProperties elas(1, 0.3);
  bie::BieElastostatic<bie::Triangle<0>, bie::Triangle<0>,bie::ElasticKernelType::H> test(elas, xy.size(1));
  std::vector<double> test_self = test.influence(source, 0, source, 0);
  std::cout << "test self effect "
             << "\n";
   for (int i = 0; i < 9; i++) {
     std::cout << test_self[i] << "\n";
   }
  ASSERT_NEAR(test_self[0], 0.656304, 1.e-3);
}

TEST(Triangle0, test_disp_1) {
    // test stress observation single element
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1.0;
    xy(2, 1) = 1.0;

    bie::Triangle<0> source;
    source.SetElement(xy);
    il::Array2D<double> xobs{1,3,1.};
    bie::Point<3> obs;
    obs.SetElement(xobs);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,bie::ElasticKernelType::T> singleT(elas, xy.size(1));
    std::cout << singleT.spatial_dimension() <<"\n";
    std::vector<double> test_displ = singleT.influence(source, 0, obs, 0);
    std::cout << "test displacement "      << "\n";
    bool test=true;
    for (int i = 0; i < 9; i++) {
        std::cout << "i: " << i << " - " << test_displ[i] << "\n";
    }
    for (int i = 8; i < 9; i++) {
        test = test && (abs(test_displ[i]+0.0206454)< 1.e-4 );
    }

    ASSERT_TRUE(test );
    std::cout <<"end test Triangle0.test_disp_1" <<"\n";
}


TEST(Triangle0, test_disp_2) {
    // test displacement observation single element - self effect
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1.0;
    xy(2, 1) = 1.0;

    bie::Triangle<0> source;
    source.SetElement(xy);
    il::Array2D<double> colots = source.collocation_points();
    il::Array2D<double> xobs{1,3,1.};
    bie::Point<3> obs;
    obs.SetElement(xobs);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,bie::ElasticKernelType::T> singleT(elas, xy.size(1));
    std::cout << singleT.spatial_dimension() <<"\n";
    std::vector<double> test_displ = singleT.influence(source, 0, obs, 0);
    std::cout << "test displacement "      << "\n";
    bool test=true;
    for (int i = 0; i < 9; i++) {
        std::cout << "i: " << i << " - " << test_displ[i] << "\n";
    }
    for (int i = 8; i < 9; i++) {
        test = test && (abs(test_displ[i]+0.0206454)< 1.e-4 );
    }

    ASSERT_TRUE(test );
    std::cout <<"end test Triangle0.test_disp_1" <<"\n";
}



TEST(Triangle0, test_stress_1) {
    // test stress observation single element
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1.0;
    xy(2, 1) = 1.0;

    bie::Triangle<0> source;
    source.SetElement(xy);
    il::Array2D<double> xobs{1,3,1.};
    bie::Point<3> obs;
    obs.SetElement(xobs);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,bie::ElasticKernelType::W> single_triangle_stress(elas, xy.size(1));
    std::vector<double> test_stress = single_triangle_stress.influence(source, 0, obs, 0);
    std::cout << "test stress "      << "\n";
    bool test=true;

    for (int i = 0; i < test_stress.size(); i++) {
        std::cout << "i: " << i << " - " << test_stress[i] << "\n";
    }
    for (int i = 16; i <  test_stress.size(); i++) {
        test = test && (abs(test_stress[i]-0.012969)< 1.e-4 );
    }

    ASSERT_TRUE(test );
    std::cout <<"end test Triangle0.test_stress_1" <<"\n";
}



TEST(Triangle0,test_Several_kernels){
    // kernel influence test for a single element
    // testing if we can have several kernels defined in the code
    il::Array2D<double> elt_coor{3,3,0.};
    elt_coor(1,0)=1.;
    elt_coor(2,1)=1.;
    bie::Triangle<0> source_elt;
    source_elt.SetElement(elt_coor);
    bie::ElasticProperties elas(1, 0.3);
    bie::BieElastostatic<bie::Triangle<0>, bie::Triangle<0>,bie::ElasticKernelType::H> testH(elas, 3);
    std::vector<double> self_effect = testH.influence(source_elt, 0, source_elt, 0);
    std::cout << "test self effect " << "\n";
     for (int i = 0; i < 9; i++) {
       std::cout << self_effect[i] << "\n";
     }

    il::Array2D<double> obs_corr{1,3,1.};
    bie::Point<3> obs_pt;
    obs_pt.SetElement(obs_corr);

    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,bie::ElasticKernelType::T> singleT(elas, 3);
    std::cout << singleT.spatial_dimension() <<"\n";
    std::vector<double> test_displ = singleT.influence(source_elt, 0, obs_pt, 0);
    std::cout << "test displacement "      << "\n";
    for (int i = 0; i < 9; i++) {
        std::cout << "i: " << i << " - " << test_displ[i] << "\n";
    }

    bie::BieElastostatic<bie::Triangle<0>, bie::Point<3>,bie::ElasticKernelType::W> singleTW(elas, 3);
    std::cout << singleTW.spatial_dimension() <<"\n";
    auto test_stres = singleTW.influence(source_elt, 0, obs_pt, 0);
    std::cout << "test stress "      << "\n";
    for (int i = 0; i < test_stres.size(); i++) {
        std::cout << test_stres[i] << "\n";
    }

    bool test_g=true;

    ASSERT_TRUE(test_g);
}