//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/math.h>

#include <src/core/BoundaryElement.h>
#include <src/elasticity/3d/BIE_elastostatic_triangle_0_impls.h>


TEST(Triangle0, test_H_1) {
    // test self-effect
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1;
    xy(2, 1) = 1.0;
    bie::Triangle<0> source;
    source.setElement(xy);
    bie::ElasticProperties elas(1, 0.3);
    bie::BIE_elastostatic<bie::Triangle<0>, bie::Triangle<0>, bie::ElasticKernelType::H> test(elas, xy.size(1));
    std::vector<double> test_self = test.influence(source, 0, source, 0);
    std::cout << "test self effect " << "\n";
    for (int i = 0; i < 9; i++) {
        std::cout << test_self[i] << "\n";
    }
    ASSERT_NEAR(test_self[0], 0.656304, 1.e-3);
}