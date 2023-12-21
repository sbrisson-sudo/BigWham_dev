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


#include "elements/point.h"

TEST(Point, test_centroid_2d) {
il::Array2D<double> xy{1, 2, 0.};
xy(0, 0) = 1.0;

bie::Point<2> pt;
pt.SetElement(xy);

auto center = pt.centroid();
std::cout << "center " << center[0] << "\n";
ASSERT_TRUE(center[0] == 1.0);
}

TEST(Point, test_collocation_2d) {
    il::Array2D<double> xy{1, 2, 0.};
    xy(0, 0) = 1.0;
    bie::Point<2> pt;
    pt.SetElement(xy);
    auto col = pt.collocation_points();
    ASSERT_TRUE(col(0,0) == 1.0 && col(0,1)==0. && col.size(0)==1);
}


TEST(Point, test_centroid_3d) {
    il::Array2D<double> xy{1, 3, 0.};
    xy(0, 0) = 1.0;
    xy(0, 2) = 3.0;
    bie::Point<3> pt;
    pt.SetElement(xy);

    auto center = pt.centroid();
    ASSERT_TRUE(center[0] == 1.0 && center[2]==3.0) ;
}


TEST(Point, test_collocation_3d) {
    il::Array2D<double> xy{1, 3, 0.};
    xy(0, 0) = 1.0;
    xy(0, 2) = 3.0;
    bie::Point<3> pt;
    pt.SetElement(xy);

    auto col = pt.collocation_points();
    ASSERT_TRUE(col(0,0) == 1.0 && col(0,2)==3.0 && col.size(0)==1) ;
}