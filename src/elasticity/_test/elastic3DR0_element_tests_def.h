//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <il/gtest/gtest/gtest.h>
#include <iostream>
#include <il/Array2D.h>
#include <src/elasticity/3d/Elastic3DR0_element.h>


//--------------------------------------------------------------------------

// Class definition
class Test3DR0Stress : public ::testing::Test{
protected:
    void TestStress(double x, double y, double  z, double  a, double  b, double G, double nu, il::Array2D<double> ReferenceStress);
};

//Method implementation
void Test3DR0Stress::TestStress(double x, double y, double  z, double  a, double  b, double G, double nu, il::Array2D<double> ReferenceStress)
{
    il::StaticArray2D<double, 3, 6> Stress;
    Stress = bie::StressesKernelR0(x, y,  z,  a,  b, G, nu);

    // stress due to displacement discontinuity DDx (shear)
    ASSERT_NEAR(ReferenceStress(0, 0), Stress(0, 0), 1.e-5);// sxx
    ASSERT_NEAR(ReferenceStress(0, 1), Stress(0, 1), 1.e-5);// syy
    ASSERT_NEAR(ReferenceStress(0, 2), Stress(0, 2), 1.e-5);// szz
    ASSERT_NEAR(ReferenceStress(0, 3), Stress(0, 3), 1.e-5);// sxy
    ASSERT_NEAR(ReferenceStress(0, 4), Stress(0, 4), 1.e-5);// sxz
    ASSERT_NEAR(ReferenceStress(0, 5), Stress(0, 5), 1.e-5);// syz

    // stress due to displacement discontinuity  DDy (shear)
    ASSERT_NEAR(ReferenceStress(1, 0), Stress(1, 0), 1.e-5);// sxx
    ASSERT_NEAR(ReferenceStress(1, 1), Stress(1, 1), 1.e-5);// syy
    ASSERT_NEAR(ReferenceStress(1, 2), Stress(1, 2), 1.e-5);// szz
    ASSERT_NEAR(ReferenceStress(1, 3), Stress(1, 3), 1.e-5);// sxy
    ASSERT_NEAR(ReferenceStress(1, 4), Stress(1, 4), 1.e-5);// sxz
    ASSERT_NEAR(ReferenceStress(1, 5), Stress(1, 5), 1.e-5);// syz

    // stress due to displacement discontinuity DDz (normal)
    ASSERT_NEAR(ReferenceStress(2, 0), Stress(2, 0), 1.e-5);// sxx
    ASSERT_NEAR(ReferenceStress(2, 1), Stress(2, 1), 1.e-5);// syy
    ASSERT_NEAR(ReferenceStress(2, 2), Stress(2, 2), 1.e-5);// szz
    ASSERT_NEAR(ReferenceStress(2, 3), Stress(2, 3), 1.e-5);// sxy
    ASSERT_NEAR(ReferenceStress(2, 4), Stress(2, 4), 1.e-5);// sxz
    ASSERT_NEAR(ReferenceStress(2, 5), Stress(2, 5), 1.e-5);// syz
}