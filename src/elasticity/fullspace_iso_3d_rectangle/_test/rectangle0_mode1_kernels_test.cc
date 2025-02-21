//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.01.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//


#include <gtest/gtest.h>

#include <il/Array2D.h>
#include <il/math.h>

#include "elements/rectangle.h"
#include "elasticity/fullspace_iso_3d_rectangle/bie_elastostatic_rectangle_0_mode1_influence.h"

TEST(Rectangle0_mode1,rect_0_H_1){
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
    bigwham::BieElastostaticModeI<bigwham::Rectangle<0>, bigwham::Rectangle<0>,bigwham::ElasticKernelType::H> test(elas, xyz.size(1));
    std::vector<double> test_self = test.influence(rec0_elt, 0, rec0_elt, 0);
    std::cout << "test self effect H "     << "\n";
    for (int i = 0; i < test_self.size(); i++) {
        std::cout << test_self[i] << "\n";
    }
    ASSERT_TRUE((rec0_elt.spatial_dimension() == 3) && (abs(test_self[0]-0.391078)<1.e-4));
}

