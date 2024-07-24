//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 17.07.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <cstdlib>
#include <memory>
#include <cmath>

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include "core/elastic_properties.h"
#include "elasticity/bie_elastostatic.h"
#include "elasticity/fullspace_iso_3d_triangle/bie_elastostatic_triangle_2_influence.h"
#include "elements/triangle.h"
#include "elements/point.h"

TEST(Triangle6, test_H_1) {
    // test self-effect triangle 0
    il::Array2D<double> xy{3, 3, 0.};
    xy(1, 0) = 1;
    xy(2, 1) = 1.0;
    bigwham::Triangle<2> source;
    source.SetElement(xy);
    bigwham::ElasticProperties elas(1, 0.3);
    bigwham::BieElastostatic<bigwham::Triangle<2>, bigwham::Triangle<2>,bigwham::ElasticKernelType::H> test(elas, xy.size(1));
    std::vector<double> test_self = test.influence(source, 0, source, 0);
    std::cout << "test self effect "     << "\n";
    for (int i = 0; i < 9; i++) {
        std::cout << test_self[i] << "\n";
    }
    ASSERT_NEAR(test_self[0], 0.635198, 1.e-3);
}
