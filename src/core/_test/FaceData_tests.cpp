//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <gtest/gtest.h>
#include <iostream>
#include <il/Array2D.h>
#include "core/oldies/FaceData.h"


//--------------------------------------------------------------------------


TEST(FaceData, centroid_1) {
    // Testing 3DR0 element 01
    il::Array2D<double> xv{il::value,{{-1.,1.,1.,-1.},{-1.,-1.,1.,1.},{0.,0.,0.,0.}}};
    il::int_t p = 0;
    bigwham::FaceData myface(xv,p);
    ASSERT_NEAR(0.,myface.collocation_points()(0,0), 1.e-5); // x
    ASSERT_NEAR(0.,myface.collocation_points()(0,1), 1.e-5); // y
    ASSERT_NEAR(0.,myface.collocation_points()(0,2), 1.e-5); // z
}

TEST(FaceData, centroid_2) {
    // Testing 3DR0 element 01
    il::Array2D<double> xv{il::value,{{-1.,1.,1.,-1.},{-1.,-1.,1.,1.},{1.,1.,1.,1.}}};
    il::int_t p = 0;
    bigwham::FaceData myface(xv,p);
    ASSERT_NEAR(0.,myface.collocation_points()(0,0), 1.e-5); // x
    ASSERT_NEAR(0.,myface.collocation_points()(0,1), 1.e-5); // y
    ASSERT_NEAR(1.,myface.collocation_points()(0,2), 1.e-5); // z
}

TEST(FaceData, centroid_3) {
    // Testing 3DR0 element 01
    il::Array2D<double> xv{il::value,{{-1.,1.,1.,-1.},{-1.,-1.,1.,1.},{1.,1.,-1.,-1.}}};
    il::int_t p = 0;
    bigwham::FaceData myface(xv,p);
    ASSERT_NEAR(0.,myface.collocation_points()(0,0), 1.e-5); // x
    ASSERT_NEAR(0.,myface.collocation_points()(0,1), 1.e-5); // y
    ASSERT_NEAR(0.,myface.collocation_points()(0,2), 1.e-5); // z
}

TEST(FaceData, centroid_4) {
    // Testing 3DR0 element 01
    il::Array2D<double> xv{il::value,{{1.,2.,2.,1.},{1.,1.,2.,2.},{1.,1.,-1.,-1.}}};
    il::int_t p = 0;
    bigwham::FaceData myface(xv,p);
    ASSERT_NEAR(1.5,myface.collocation_points()(0,0), 1.e-5); // x
    ASSERT_NEAR(1.5,myface.collocation_points()(0,1), 1.e-5); // y
    ASSERT_NEAR(0.,myface.collocation_points()(0,2), 1.e-5);  // z
}
