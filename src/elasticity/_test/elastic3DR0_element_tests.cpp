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


TEST(R0, ip_expressions) {
  il::Array2D<double> ReferenceValues{19, 1, 0.};
  double x, y, z, xi, eta, res;
  x = 1.12;
  y = 1.3;
  z = 2.14;
  xi = 1.13;
  eta = 2.21;

  // Reference values from Mathematica with 18 decimals
  ReferenceValues(0,0) = 0.3474600029108964;        //ip1
  ReferenceValues(1,0) = 0.8396116951656556;        //ip2
  ReferenceValues(2,0) = -0.0018285920369267496;    //ip3
  ReferenceValues(3,0) = -0.003038013695351292;     //ip11
  ReferenceValues(4,0) = 0.4300210242418253;        //ip12
  ReferenceValues(5,0) = 0.6501349308051904;        //ip13
  ReferenceValues(6,0) = -0.1690021989608522;       //ip22
  ReferenceValues(7,0) = 0.3974337426112349;        //ip23
  ReferenceValues(8,0) = 0.0015780977665905236;     //ip33
  ReferenceValues(9,0) = 0.30378652217128777;       //ip111
  ReferenceValues(10,0) = 0.0007951866271715226;    //ip112
  ReferenceValues(11,0) = 0.0031773358634533603;    //ip113
  ReferenceValues(12,0) = 0.07236198307261016;      //ip122
  ReferenceValues(13,0) = -0.1701699382147095;      //ip123
  ReferenceValues(14,0) = -0.37614850524389787;     //ip133
  ReferenceValues(15,0) = 0.12871603714845617;      //ip222
  ReferenceValues(16,0) = 0.1340455199049876;       //ip223
  ReferenceValues(17,0) = -0.12951122377562768;     //ip233
  ReferenceValues(18,0) = -0.0019957694764766985;   //ip333

  ASSERT_NEAR(ReferenceValues(0,0), bie::ip1(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(1,0), bie::ip2(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(2,0), bie::ip3(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(3,0), bie::ip11(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(4,0), bie::ip12(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(5,0), bie::ip13(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(6,0), bie::ip22(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(7,0), bie::ip23(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(8,0), bie::ip33(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(9,0), bie::ip111(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(10,0), bie::ip112(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(11,0), bie::ip113(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(12,0), bie::ip122(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(13,0), bie::ip123(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(14,0), bie::ip133(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(15,0), bie::ip222(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(16,0), bie::ip223(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(17,0), bie::ip233(x, y, z, xi, eta), 1.e-5);
  ASSERT_NEAR(ReferenceValues(18,0), bie::ip333(x, y, z, xi, eta), 1.e-5);
}

TEST(R0, stress_expressions) {
    il::StaticArray2D<double, 3, 6> Stress;
    il::Array2D<double> ReferenceStress{3, 6, 0.};
    double x, y, z, a, b, G, nu;
    x = 1.12;
    y = 1.3;
    z = 2.14;
    a = 1.;
    b = 2.;
    G = 10.;
    nu = 0.3;
    Stress = bie::StressesKernelR0(x, y,  z,  a,  b, G, nu);

    // Reference values from Mathematica with 18 decimals
    // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

    // stress due to displacement discontinuity DDx (shear)
    ReferenceStress(0, 0) = 0.025135933881063915;  // sxx
    ReferenceStress(0, 1) = 0.07714206469390152;   // syy
    ReferenceStress(0, 2) = 0.807226930217386;     // szz
    ReferenceStress(0, 3) = 0.06764260339990966;   // sxy
    ReferenceStress(0, 4) = 0.1182294751296132;    // sxz
    ReferenceStress(0, 5) = 0.220426565874335;     // syz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceStress(1, 0) = 0.04428228498873277;    // sxx
    ReferenceStress(1, 1) = 0.1710180878105341;     // syy
    ReferenceStress(1, 2) = 0.3920679058913325;     // szz
    ReferenceStress(1, 3) = 0.1121230234936069;     // sxy
    ReferenceStress(1, 4) = 0.2204265658743351;     // sxz
    ReferenceStress(1, 5) = 0.05116011738468318;    // syz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceStress(2, 0) = 0.0975375342845718;     // sxx
    ReferenceStress(2, 1) = 0.02894407986186851;    // syy
    ReferenceStress(2, 2) = 0.989125823417809;      // szz
    ReferenceStress(2, 3) = 0.2102150629008792;     // sxy
    ReferenceStress(2, 4) = 0.807226930217385;      // sxz
    ReferenceStress(2, 5) = 0.3920679058913323;     // syz

    ASSERT_NEAR(ReferenceStress(0, 0), Stress(0, 0), 1.e-5);
    ASSERT_NEAR(ReferenceStress(0, 1), Stress(0, 1), 1.e-5);
    ASSERT_NEAR(ReferenceStress(0, 2), Stress(0, 2), 1.e-5);
    ASSERT_NEAR(ReferenceStress(0, 3), Stress(0, 3), 1.e-5);
    ASSERT_NEAR(ReferenceStress(0, 4), Stress(0, 4), 1.e-5);
    ASSERT_NEAR(ReferenceStress(0, 5), Stress(0, 5), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 0), Stress(1, 0), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 1), Stress(1, 1), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 2), Stress(1, 2), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 3), Stress(1, 3), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 4), Stress(1, 4), 1.e-5);
    ASSERT_NEAR(ReferenceStress(1, 5), Stress(1, 5), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 0), Stress(2, 0), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 1), Stress(2, 1), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 2), Stress(2, 2), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 3), Stress(2, 3), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 4), Stress(2, 4), 1.e-5);
    ASSERT_NEAR(ReferenceStress(2, 5), Stress(2, 5), 1.e-5);
}

TEST(R0, displacement_expressions) {
    il::Array2D<double> Displacement;
    il::Array2D<double> ReferenceDisplacement{3, 3, 0.};
    double x, y, z, a, b, nu;
    x = 1.12;
    y = 1.3;
    z = 2.14;
    a = 1.;
    b = 2.;
    nu = 0.3;
    Displacement = bie::DisplacementKernelR0( x, y, z, a, b, nu);
    // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
    //   0      -> |       Ux,            Ux,             Ux            |
    //   1      -> |       Uy,            Uy,             Uy            |
    //   2      -> |       Uz,            Uz,             Uz            |

    // Reference values from Mathematica with 18 decimals
    // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

    // stress due to displacement discontinuity DDx (shear)
    ReferenceDisplacement(0, 0) = -0.03892389878412812;    // ux
    ReferenceDisplacement(1, 0) = -0.01092630818159775;    // uy
    ReferenceDisplacement(2, 0) = -0.0456133521126238;     // uz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceDisplacement(0, 1) = -0.01092630818159775;     // ux
    ReferenceDisplacement(1, 1) = -0.04055468222934548;     // uy
    ReferenceDisplacement(2, 1) = -0.03183576694099982;     // uz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceDisplacement(0, 2) = -0.02924589971874653;     // ux
    ReferenceDisplacement(1, 2) = -0.01815531445891874;     // uy
    ReferenceDisplacement(2, 2) = -0.1086065957871454;      // uz


    ASSERT_NEAR(ReferenceDisplacement(0, 0), Displacement(0, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(0, 1), Displacement(0, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(0, 2), Displacement(0, 2), 1.e-5);

    ASSERT_NEAR(ReferenceDisplacement(1, 0), Displacement(1, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(1, 1), Displacement(1, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(1, 2), Displacement(1, 2), 1.e-5);

    ASSERT_NEAR(ReferenceDisplacement(2, 0), Displacement(2, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(2, 1), Displacement(2, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(2, 2), Displacement(2, 2), 1.e-5);

}
