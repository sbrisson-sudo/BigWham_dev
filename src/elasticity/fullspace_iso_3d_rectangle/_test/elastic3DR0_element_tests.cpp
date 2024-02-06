//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <gtest/gtest.h>
#include <iostream>
#include <il/Array2D.h>
#include "elasticity/fullspace_iso_3d_rectangle//elastic_3dR0_common.h"
#include "elasticity/fullspace_iso_3d_rectangle/elastic_3dR0_element.h"
#include "elastic3DR0_element_tests_def.h"

//--------------------------------------------------------------------------


TEST(R0, ip_expressions) {
  il::Array2D<double> ReferenceValues{19, 1, 0.};
  double x, y, z, xi, eta;
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


TEST(R0, Ip_special_expressions) {
    il::Array2D<double> ReferenceValues{2, 1, 0.};
    double x, y, a, b;
    a = 2.000;
    b = 3.000;

    // Reference values from Mathematica with 18 decimals
    ReferenceValues(0,0) = -0.14308728328344156;       //ip33 limits z->0 and x->a
    ReferenceValues(1,0) = -0.0921292729553324;        //ip33 limits z->0 and y->b

    ASSERT_NEAR(ReferenceValues(0,0), bie::Ip33_lim_z_to_0_and_x_to_a (x=a, y=6.000, a, b), 1.e-5);
    ASSERT_NEAR(ReferenceValues(1,0), bie::Ip33_lim_z_to_0_and_y_to_b (x=6.000, y=b, a, b), 1.e-5);
}

TEST(R0, ip_cheking_limits) {
    il::Array2D<double> ReferenceValues{3, 1, 0.};
    double x, y, z=0.000, xi, eta;
    x = 4.000;
    y = 6.000;
    xi = 2.000;
    eta = 3.000;

    // Reference values from Mathematica with 18 decimals

    ReferenceValues(0,0) = 0.0839748528310782;       //ip11, z->0
    ReferenceValues(1,0) = 0.1484332679249236;       //ip22, z->0
    ReferenceValues(2,0) = 0.6009252125773314;       //ip33, z->0

    ASSERT_NEAR(ReferenceValues(0,0), bie::ip11(x, y, z, xi, eta), 1.e-5);
    ASSERT_NEAR(ReferenceValues(1,0), bie::ip22(x, y, z, xi, eta), 1.e-5);
    ASSERT_NEAR(ReferenceValues(2,0), bie::ip33(x, y, z, xi, eta), 1.e-5);
}

TEST_F(Test3DR0Stress, stress_expressions_at_non_zero_z) {
    double x, y, z, a, b, G, nu;
    x = 1.12;
    y = 1.3;
    z = 2.14;
    a = 1.;
    b = 2.;
    G = 10.;
    nu = 0.3;

    il::Array2D<double> ReferenceStress{3,6,0.};
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

    TestStress(x, y, z, a, b, G, nu, ReferenceStress);
}


TEST_F(Test3DR0Stress, stress_expressions_at_zero_z) {
    double x, y, z, a, b, G, nu;
    x = 0.9;
    y = 3.1000000;
    z = 0.0000000;
    a = 1.;
    b = 2.;
    G = 200.;
    nu = 0.3;

    il::Array2D<double> ReferenceStress{3,6,0.};
    // Reference values from Mathematica with 18 decimals
    // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

    ReferenceStress(0, 0) = 0.0;                    // sxx
    ReferenceStress(0, 1) = 0.0;                    // syy
    ReferenceStress(0, 2) = 0.0;                    // szz
    ReferenceStress(0, 3) = 0.0;                    // sxy
    ReferenceStress(0, 4) = -6.572822442880559;     // sxz
    ReferenceStress(0, 5) = -2.98463799233973;      // syz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceStress(1, 0) = 0.0;                     // sxx
    ReferenceStress(1, 1) = 0.0;                     // syy
    ReferenceStress(1, 2) = 0.0;                     // szz
    ReferenceStress(1, 3) = 0.0;                     // sxy
    ReferenceStress(1, 4) = -2.98463799233973;       // sxz
    ReferenceStress(1, 5) = -13.83986218452465;      // syz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceStress(2, 0) = -4.761276075329918;      // sxx
    ReferenceStress(2, 1) = -14.45066239752203;     // syy
    ReferenceStress(2, 2) = -12.007461545532491;     // szz
    ReferenceStress(2, 3) = -3.9795173231196403;      // sxy
    ReferenceStress(2, 4) = 0.0;                     // sxz
    ReferenceStress(2, 5) = 0.0;                     // syz

    TestStress(x, y, z, a, b, G, nu, ReferenceStress);
}

TEST_F(Test3DR0Stress, stress_expressions_at_zero_z_and_x_to_a) {
    double x, y, z, a, b, G, nu;
    x = 1.;
    y = 3.1000000;
    z = 0.0000000;
    a = 1.;
    b = 2.;
    G = 200.;
    nu = 0.3;

    il::Array2D<double> ReferenceStress{3,6,0.};
    // Reference values from Mathematica with 18 decimals
    // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

    ReferenceStress(0, 0) = 0.0;                    // sxx
    ReferenceStress(0, 1) = 0.0;                    // syy
    ReferenceStress(0, 2) = 0.0;                    // szz
    ReferenceStress(0, 3) = 0.0;                    // sxy
    ReferenceStress(0, 4) = -6.433376219207982;     // sxz
    ReferenceStress(0, 5) = -3.120220980453658;      // syz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceStress(1, 0) = 0.0;                     // sxx
    ReferenceStress(1, 1) = 0.0;                     // syy
    ReferenceStress(1, 2) = 0.0;                     // szz
    ReferenceStress(1, 3) = 0.0;                     // sxy
    ReferenceStress(1, 4) = -3.120220980453658;       // sxz
    ReferenceStress(1, 5) = -12.9098547996334;      // syz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceStress(2, 0) = -4.785044563092725;      // sxx
    ReferenceStress(2, 1) = -13.42034933699328;     // syy
    ReferenceStress(2, 2) = -11.37837118755376;     // szz
    ReferenceStress(2, 3) = -4.160294640604879;      // sxy
    ReferenceStress(2, 4) = 0.0;                     // sxz
    ReferenceStress(2, 5) = 0.0;                     // syz

    TestStress(x, y, z, a, b, G, nu, ReferenceStress);
}


TEST_F(Test3DR0Stress, stress_expressions_at_zero_z_and_y_to_b) {
    double x, y, z, a, b, G, nu;
    x = 1.50000;
    y = 2.00000;
    z = 0.0000000;
    a = 1.00000;
    b = 2.00000;
    G = 200.;
    nu = 0.3;

    il::Array2D<double> ReferenceStress{3,6,0.};
    // Reference values from Mathematica with 18 decimals
    // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

    ReferenceStress(0, 0) = 0.0;                    // sxx
    ReferenceStress(0, 1) = 0.0;                    // syy
    ReferenceStress(0, 2) = 0.0;                    // szz
    ReferenceStress(0, 3) = 0.0;                    // sxy
    ReferenceStress(0, 4) = -35.79423536018405;     // sxz
    ReferenceStress(0, 5) = -10.66745173504686;      // syz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceStress(1, 0) = 0.0;                     // sxx
    ReferenceStress(1, 1) = 0.0;                     // syy
    ReferenceStress(1, 2) = 0.0;                     // szz
    ReferenceStress(1, 3) = 0.0;                     // sxy
    ReferenceStress(1, 4) = -10.66745173504686;       // sxz
    ReferenceStress(1, 5) = -23.87911771266086;      // syz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceStress(2, 0) = -36.02498968164837;      // sxx
    ReferenceStress(2, 1) = -20.13816615161743;     // syy
    ReferenceStress(2, 2) = -35.10197239579112;     // szz
    ReferenceStress(2, 3) = -14.223268980062478;      // sxy
    ReferenceStress(2, 4) = 0.0;                     // sxz
    ReferenceStress(2, 5) = 0.0;                     // syz

    TestStress(x, y, z, a, b, G, nu, ReferenceStress);
}

TEST(R0, singular_stress){
    double x, y, z, a, b;
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 1., y = 2., z = 0., a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -1., y = 2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -1., y = -2.,z = 0., a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 1., y = -2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 1., y = 1., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -1., y = 1., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -1., y = -1., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 1., y = -1., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 1., y = 10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -1., y = 10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -1., y = -10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 1., y = -10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 0.5, y = 2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -0.5, y = 2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = -0.5, y = -2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 0.5, y = -2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 10., y = 2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -10., y = 2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -10., y = -2., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 10., y = -2.,z = 0., a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 0., y = 0.,z = 0., a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 10., y = 10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -10., y = 10.,z = 0., a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = -10., y = -10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 10., y = -10., z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(false, bie::is_stress_singular_at_given_location(x = 0.9999, y = -1.99999, z = 0.,a = 1., b = 2., false));
    ASSERT_EQ(true, bie::is_stress_singular_at_given_location(x = 0.999999999999999999999, y = -1.999999999999999999999,z = 0., a = 1., b = 2., false));
};

TEST(R0, displ_expr) {
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

    x = 0.00;
    y = 0.00;
    z = 1.00;
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
    ReferenceDisplacement(0, 0) = -0.1251318438095986;    // ux
    ReferenceDisplacement(1, 0) = 0.;    // uy
    ReferenceDisplacement(2, 0) = 0.;     // uz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceDisplacement(0, 1) = 0.;     // ux
    ReferenceDisplacement(1, 1) = -0.180824472469147;     // uy
    ReferenceDisplacement(2, 1) = 0.;     // uz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceDisplacement(0, 2) = 0.;     // ux
    ReferenceDisplacement(1, 2) = 0.;     // uy
    ReferenceDisplacement(2, 2) = -0.347902358447792;      // uz


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

TEST(R0, displ_expr_on_the_elem_plane) {
    il::Array2D<double> Displacement;
    il::Array2D<double> ReferenceDisplacement{3, 3, 0.};
    double x, y, z, a, b, nu;
    x = 1.5;
    y = 2.5;
    z = 0.;
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
    ReferenceDisplacement(0, 0) = 0.;                       // ux
    ReferenceDisplacement(1, 0) = 0.;                       // uy
    ReferenceDisplacement(2, 0) = -0.01956059198093905;     // uz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceDisplacement(0, 1) = 0.;                       // ux
    ReferenceDisplacement(1, 1) = 0.;                       // uy
    ReferenceDisplacement(2, 1) = -0.02300029715436082;     // uz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceDisplacement(0, 2) = 0.019560591980939054;     // ux
    ReferenceDisplacement(1, 2) = 0.02300029715436082;      // uy
    ReferenceDisplacement(2, 2) = 0.;                       // uz


    ASSERT_NEAR(ReferenceDisplacement(0, 0), Displacement(0, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(0, 1), Displacement(0, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(0, 2), Displacement(0, 2), 1.e-5);

    ASSERT_NEAR(ReferenceDisplacement(1, 0), Displacement(1, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(1, 1), Displacement(1, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(1, 2), Displacement(1, 2), 1.e-5);

    ASSERT_NEAR(ReferenceDisplacement(2, 0), Displacement(2, 0), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(2, 1), Displacement(2, 1), 1.e-5);
    ASSERT_NEAR(ReferenceDisplacement(2, 2), Displacement(2, 2), 1.e-5);

    x = 0.25;
    y = 1.10;
    z = 0.00;
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
    ReferenceDisplacement(0, 0) = -0.5;    // ux
    ReferenceDisplacement(1, 0) = 0.;    // uy
    ReferenceDisplacement(2, 0) = -0.01895705258167849;     // uz

    // stress due to displacement discontinuity  DDy (shear)
    ReferenceDisplacement(0, 1) = 0.;     // ux
    ReferenceDisplacement(1, 1) = -0.5;     // uy
    ReferenceDisplacement(2, 1) = -0.028587157479065407;     // uz

    // stress due to displacement discontinuity DDz (normal)
    ReferenceDisplacement(0, 2) = 0.01895705258167849;     // ux
    ReferenceDisplacement(1, 2) = 0.0285871574790654;     // uy
    ReferenceDisplacement(2, 2) = -0.5;      // uz


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
