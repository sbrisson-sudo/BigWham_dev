//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.01.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//



#include <gtest/gtest.h>

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include "elasticity/fullspace_iso_3d_rectangle/bie_elastostatic_rectangle_0_influence.h"
#include "elements/rectangle.h"

TEST(Rectangle0,rect_0_H_1){
    double hx = 1.;
    double hy = 2.;
    il::Array2D<double> xyz{4, 3, 0.0};
    xyz(1, 0) = hx;
    xyz(2, 0) = hx;
    xyz(2, 1) = hy;
    xyz(3, 1) = hy;
    bigwham::Rectangle<0> rec0_elt;
    rec0_elt.SetElement(xyz);
//test self effect H kernel
    bigwham::ElasticProperties elas(1, 0.3);
    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Rectangle<0>,bigwham::ElasticKernelType::H> test(elas, xyz.size(1));
    std::vector<double> test_self = test.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test self effect H "     << "\n";
    for (int i = 0; i < test_self.size(); i++) {
        std::cout << test_self[i] << "\n";
    }

    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) && (abs(test_self[8]-0.391078)<1.e-4));
}


TEST(Rectangle0,rect_0_T_1){
    double hx = 1.;
    double hy = 2.;
    il::Array2D<double> xyz{4, 3, 0.0};
    xyz(1, 0) = hx;
    xyz(2, 0) = hx;
    xyz(2, 1) = hy;
    xyz(3, 1) = hy;
    bigwham::Rectangle<0> rec0_elt;
    rec0_elt.SetElement(xyz);
//test self effect T kernel
    bigwham::ElasticProperties elas(1, 0.3);
    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Rectangle<0>,bigwham::ElasticKernelType::T> test(elas, xyz.size(1));
    std::vector<double> test_self = test.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test self effect T"     << "\n";
    for (int i = 0; i < test_self.size(); i++) {
        std::cout << test_self[i] << "\n";
    }

    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) && (test_self[0]==-0.5) && (test_self[4]==-0.5) && (test_self[8]==-0.5) );
}


TEST(Rectangle0,rect_0_disp_obs){
    double hx = 1.;
    double hy = 2.;
    il::Array2D<double> xyz{4, 3, 0.0};
    xyz(1, 0) = hx;
    xyz(2, 0) = hx;
    xyz(2, 1) = hy;
    xyz(3, 1) = hy;
    bigwham::Rectangle<0> rec0_elt;
    rec0_elt.SetElement(xyz);
    bigwham::ElasticProperties elas(1, 0.3);
    il::Array2D<double> xobs{1,3,1.};
    bigwham::Point<3> obs;
    obs.SetElement(xobs);
    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Point<3>,bigwham::ElasticKernelType::T> test_disp(elas, xyz.size(1));
    std::vector<double> test_disp_ = test_disp.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test displacement observation - "  <<  test_disp_.size()  << "\n";
    for (int i = 0; i < test_disp_.size(); i++) {
        std::cout << "i : " << test_disp_[i] << "\n";
    }

    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) );
}


TEST(Rectangle0,rect_0_stress_obs){
    double hx = 1.;
    double hy = 2.;
    il::Array2D<double> xyz{4, 3, 0.0};
    xyz(1, 0) = hx;
    xyz(2, 0) = hx;
    xyz(2, 1) = hy;
    xyz(3, 1) = hy;
    bigwham::Rectangle<0> rec0_elt;
    rec0_elt.SetElement(xyz);
    bigwham::ElasticProperties elas(1, 0.3);
    il::Array2D<double> xobs{1,3,1.};
    bigwham::Point<3> obs;
    obs.SetElement(xobs);
    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Point<3>,bigwham::ElasticKernelType::W> test_stress(elas, xyz.size(1));
    std::vector<double> test_stress_ = test_stress.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test stress - "  <<  test_stress_.size()  << "\n";
    for (int i = 0; i < test_stress_.size(); i++) {
        std::cout << "i : " << test_stress_[i] << "\n";
    }

    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) );
}

// with 2 kernels defined

TEST(Rectangle0,rect_0_disp_and_stress_obs){
    double hx = 1.;
    double hy = 2.;
    il::Array2D<double> xyz{4, 3, 0.0};
    xyz(1, 0) = hx;
    xyz(2, 0) = hx;
    xyz(2, 1) = hy;
    xyz(3, 1) = hy;
    bigwham::Rectangle<0> rec0_elt;
    rec0_elt.SetElement(xyz);
    bigwham::ElasticProperties elas(1, 0.3);
    il::Array2D<double> xobs{1,3,1.};
    bigwham::Point<3> obs;
    obs.SetElement(xobs);

    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Point<3>,bigwham::ElasticKernelType::W> test_stress(elas, xyz.size(1));
    std::vector<double> test_stress_ = test_stress.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test stress - "  <<  test_stress_.size()  << "\n";
    for (int i = 0; i < test_stress_.size(); i++) {
        std::cout << "i : " << test_stress_[i] << "\n";
    }

    bigwham::BieElastostatic<bigwham::Rectangle<0>, bigwham::Point<3>,bigwham::ElasticKernelType::T> test_disp(elas, xyz.size(1));
    std::vector<double> test_disp_ = test_disp.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test displacement observation - "  <<  test_disp_.size()  << "\n";
    for (int i = 0; i < test_disp_.size(); i++) {
        std::cout << "i : " << test_disp_[i] << "\n";
    }


    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) );
}


