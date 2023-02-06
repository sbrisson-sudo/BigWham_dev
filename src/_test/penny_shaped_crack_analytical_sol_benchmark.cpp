//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include "_test/penny_shaped_crack_analytical_sol.h"

// ---------------------------  part 1 -------------------------------------
/*
 * *****************
 * TEST FOR Z > 0  *
 * *****************
 *
 * DISPLACEMENTS
 * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
 * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
 * expressions valid for any z
 *
 */

// normal load

TEST(radial_analytical_positive_Z, U_nl_){
    double x = 1.2, y = 1.4, z = 1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double pz = 160.;
    double numerical_Ux_nl_ = bie::Ux_nl_(x, y, z, a, G, nu, pz);
    double analytical_Ux_nl_ = 3.785679043102655;
    ASSERT_NEAR(numerical_Ux_nl_, analytical_Ux_nl_, 1.e-5);

    double numerical_Uy_nl_ = bie::Uy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Uy_nl_ = 4.416625550286431;
    ASSERT_NEAR(numerical_Uy_nl_, analytical_Uy_nl_, 1.e-5);

    double numerical_Uz_nl_ = bie::Uz_nl_(x, y, z, a, G, nu, pz);
    double analytical_Uz_nl_ = 11.38838742452168;
    ASSERT_NEAR(numerical_Uz_nl_, analytical_Uz_nl_, 1.e-5);
}
// shear load

TEST(radial_analytical_positive_Z, U_sl_){
    double x = 1.2, y = 1.4, z = 1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double px = 155., py = 230.;
    double numerical_Ux_sl_ = bie::Ux_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Ux_sl_ = 14.41599475709211;
    ASSERT_NEAR(numerical_Ux_sl_, analytical_Ux_sl_, 1.e-5);

    double numerical_Uy_sl_ = bie::Uy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Uy_sl_ = 18.0412232699112;
    ASSERT_NEAR(numerical_Uy_sl_, analytical_Uy_sl_, 1.e-5);

    double numerical_Uz_sl_ = bie::Uz_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Uz_sl_ = 22.42071601458012;
    ASSERT_NEAR(numerical_Uz_sl_, analytical_Uz_sl_, 1.e-5);
}

/*
 * *****************
 * TEST FOR Z > 0  *
 * *****************
 *
 * STRESSES
 * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
 * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
 * expressions valid for any z
 *
 */

// normal load

TEST(radial_analytical_positive_Z, Stress_nl_){
    double x = 1.2, y = 1.4, z = 1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double pz = 160.;
    double numerical_Sig_xx_nl_ = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_xx_nl_ = -2.0657867576932993;
    ASSERT_NEAR(numerical_Sig_xx_nl_, analytical_Sig_xx_nl_, 1.e-5);

    double numerical_Sig_yy_nl_ = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_yy_nl_ = -3.656171942425667;
    ASSERT_NEAR(numerical_Sig_yy_nl_, analytical_Sig_yy_nl_, 1.e-5);

    double numerical_Sig_zz_nl_ = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zz_nl_ = 0.543869725910124;
    ASSERT_NEAR(numerical_Sig_zz_nl_, analytical_Sig_zz_nl_, 1.e-5);

    double numerical_Sig_xy_nl_ = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_xy_nl_ = -5.13816751990457;
    ASSERT_NEAR(numerical_Sig_xy_nl_, analytical_Sig_xy_nl_, 1.e-5);

    double numerical_Sig_zx_nl_ = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zx_nl_ = -5.654491186524116;
    ASSERT_NEAR(numerical_Sig_zx_nl_, analytical_Sig_zx_nl_, 1.e-5);

    double numerical_Sig_zy_nl_ = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zy_nl_ = -6.596906384278137;
    ASSERT_NEAR(numerical_Sig_zy_nl_, analytical_Sig_zy_nl_, 1.e-5);

}

// shear load



TEST(radial_analytical_positive_Z, Stress_sl_){
    double x = 1.2, y = 1.4, z = 1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double px = 155., py = 230.;
    double numerical_Sig_xx_sl_ = bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_xx_sl_ = -6.376791154478333;
    ASSERT_NEAR(numerical_Sig_xx_sl_, analytical_Sig_xx_sl_, 1.e-5);

    double numerical_Sig_yy_sl_ = bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_yy_sl_ = -9.2234549711472;
    ASSERT_NEAR(numerical_Sig_yy_sl_, analytical_Sig_yy_sl_, 1.e-5);

    double numerical_Sig_zz_sl_ = bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zz_sl_ = -17.600989722758893;
    ASSERT_NEAR(numerical_Sig_zz_sl_, analytical_Sig_zz_sl_, 1.e-5);

    double numerical_Sig_xy_sl_ = bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_xy_sl_ = -9.89471480535422;
    ASSERT_NEAR(numerical_Sig_xy_sl_, analytical_Sig_xy_sl_, 1.e-5);

    double numerical_Sig_zx_sl_ = bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zx_sl_ = -11.60516358682423;
    ASSERT_NEAR(numerical_Sig_zx_sl_, analytical_Sig_zx_sl_, 1.e-5);

    double numerical_Sig_zy_sl_ = bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zy_sl_ = -12.62186684505416;
    ASSERT_NEAR(numerical_Sig_zy_sl_, analytical_Sig_zy_sl_, 1.e-5);
}

// ---------------------------  part 2 -------------------------------------

/*
 * *****************
 * TEST FOR Z < 0  *
 * *****************
 *
 * DISPLACEMENTS
 * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
 * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
 * expressions valid for any z
 *
 */

// normal load

TEST(radial_analytical_negative_Z, U_nl_){
    double x = 1.2, y = 1.4, z = -1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double pz = 160.;

    double numerical_Ux_nl_ = bie::Ux_nl_(x, y, z, a, G, nu, pz);
    double analytical_Ux_nl_ = +3.7856790431026552;
    ASSERT_NEAR(numerical_Ux_nl_, analytical_Ux_nl_, 1.e-5);

    double numerical_Uy_nl_ = bie::Uy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Uy_nl_ = +4.416625550286431;
    ASSERT_NEAR(numerical_Uy_nl_, analytical_Uy_nl_, 1.e-5);

    double numerical_Uz_nl_ = bie::Uz_nl_(x, y, z, a, G, nu, pz);
    double analytical_Uz_nl_ = -11.38838742452168;
    ASSERT_NEAR(numerical_Uz_nl_, analytical_Uz_nl_, 1.e-5);

    x = -1.2; y = 1.4;  z = -1.5;
    numerical_Ux_nl_ = bie::Ux_nl_(x, y, z, a, G, nu, pz);
    analytical_Ux_nl_ = -3.7856790431026552;
    ASSERT_NEAR(numerical_Ux_nl_, analytical_Ux_nl_, 1.e-5);

    numerical_Uy_nl_ = bie::Uy_nl_(x, y, z, a, G, nu, pz);
    analytical_Uy_nl_ = +4.416625550286431;
    ASSERT_NEAR(numerical_Uy_nl_, analytical_Uy_nl_, 1.e-5);

    numerical_Uz_nl_ = bie::Uz_nl_(x, y, z, a, G, nu, pz);
    analytical_Uz_nl_ = -11.38838742452168;
    ASSERT_NEAR(numerical_Uz_nl_, analytical_Uz_nl_, 1.e-5);

    x = -1.2; y = -1.4; z = -1.5;
    numerical_Ux_nl_ = bie::Ux_nl_(x, y, z, a, G, nu, pz);
    analytical_Ux_nl_ = -3.7856790431026552;
    ASSERT_NEAR(numerical_Ux_nl_, analytical_Ux_nl_, 1.e-5);

    numerical_Uy_nl_ = bie::Uy_nl_(x, y, z, a, G, nu, pz);
    analytical_Uy_nl_ = -4.416625550286431;
    ASSERT_NEAR(numerical_Uy_nl_, analytical_Uy_nl_, 1.e-5);

    numerical_Uz_nl_ = bie::Uz_nl_(x, y, z, a, G, nu, pz);
    analytical_Uz_nl_ = -11.38838742452168;
    ASSERT_NEAR(numerical_Uz_nl_, analytical_Uz_nl_, 1.e-5);

    x =  1.2; y = -1.4; z = -1.5;
    numerical_Ux_nl_ = bie::Ux_nl_(x, y, z, a, G, nu, pz);
    analytical_Ux_nl_ = 3.7856790431026552;
    ASSERT_NEAR(numerical_Ux_nl_, analytical_Ux_nl_, 1.e-5);

    numerical_Uy_nl_ = bie::Uy_nl_(x, y, z, a, G, nu, pz);
    analytical_Uy_nl_ = -4.416625550286431;
    ASSERT_NEAR(numerical_Uy_nl_, analytical_Uy_nl_, 1.e-5);

    numerical_Uz_nl_ = bie::Uz_nl_(x, y, z, a, G, nu, pz);
    analytical_Uz_nl_ = -11.38838742452168;
    ASSERT_NEAR(numerical_Uz_nl_, analytical_Uz_nl_, 1.e-5);

}

// shear load

TEST(radial_analytical_negative_Z, U_sl_){
    double x = 1.2, y = 1.4, z = -1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double px = 155., py = 230.;

    double numerical_Ux_sl_ = bie::Ux_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Ux_sl_ = -14.41599475709211;
    ASSERT_NEAR(numerical_Ux_sl_, analytical_Ux_sl_, 1.e-5);

    double numerical_Uy_sl_ = bie::Uy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Uy_sl_ = -18.0412232699112;
    ASSERT_NEAR(numerical_Uy_sl_, analytical_Uy_sl_, 1.e-5);

    double numerical_Uz_sl_ = bie::Uz_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Uz_sl_ = 22.42071601458012;
    ASSERT_NEAR(numerical_Uz_sl_, analytical_Uz_sl_, 1.e-5);

    x =  -1.2; y = 1.4; z = -1.5;
    numerical_Ux_sl_ = bie::Ux_sl_(x, y, z, a, G, nu, px, py);
    analytical_Ux_sl_ = -1.0266085318973426;
    ASSERT_NEAR(numerical_Ux_sl_, analytical_Ux_sl_, 1.e-5);

    numerical_Uy_sl_ = bie::Uy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uy_sl_ = -9.01794124858422;
    ASSERT_NEAR(numerical_Uy_sl_, analytical_Uy_sl_, 1.e-5);

    numerical_Uz_sl_ = bie::Uz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uz_sl_ = 6.002396413352189;
    ASSERT_NEAR(numerical_Uz_sl_, analytical_Uz_sl_, 1.e-5);

    x =  -1.2; y = -1.4; z = -1.5;
    numerical_Ux_sl_ = bie::Ux_sl_(x, y, z, a, G, nu, px, py);
    analytical_Ux_sl_ = -14.41599475709211;
    ASSERT_NEAR(numerical_Ux_sl_, analytical_Ux_sl_, 1.e-5);

    numerical_Uy_sl_ = bie::Uy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uy_sl_ = -18.0412232699112;
    ASSERT_NEAR(numerical_Uy_sl_, analytical_Uy_sl_, 1.e-5);

    numerical_Uz_sl_ = bie::Uz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uz_sl_ = -22.42071601458012;
    ASSERT_NEAR(numerical_Uz_sl_, analytical_Uz_sl_, 1.e-5);

    x =  1.2; y = -1.4; z = -1.5;
    numerical_Ux_sl_ = bie::Ux_sl_(x, y, z, a, G, nu, px, py);
    analytical_Ux_sl_ = -1.026608531897343;
    ASSERT_NEAR(numerical_Ux_sl_, analytical_Ux_sl_, 1.e-5);

    numerical_Uy_sl_ = bie::Uy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uy_sl_ = -9.017941248584219;
    ASSERT_NEAR(numerical_Uy_sl_, analytical_Uy_sl_, 1.e-5);

    numerical_Uz_sl_ = bie::Uz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Uz_sl_ = -6.002396413352189;
    ASSERT_NEAR(numerical_Uz_sl_, analytical_Uz_sl_, 1.e-5);

}

/*
 * *****************
 * TEST FOR Z < 0  *
 * *****************
 *
 * STRESSES
 * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
 * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
 * expressions valid for any z
 *
 */

 //normal load

TEST(radial_analytical_negative_Z, Stress_nl_){
    double x = 1.2, y = 1.4, z = -1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double pz = 160.;

    double numerical_Sig_xx_nl_ = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_xx_nl_ = -2.0657867576932993;
    ASSERT_NEAR(numerical_Sig_xx_nl_, analytical_Sig_xx_nl_, 1.e-5);

    double numerical_Sig_yy_nl_ = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_yy_nl_ = -3.656171942425667;
    ASSERT_NEAR(numerical_Sig_yy_nl_, analytical_Sig_yy_nl_, 1.e-5);

    double numerical_Sig_zz_nl_ = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zz_nl_ = 0.543869725910124;
    ASSERT_NEAR(numerical_Sig_zz_nl_, analytical_Sig_zz_nl_, 1.e-5);

    double numerical_Sig_xy_nl_ = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_xy_nl_ = -5.13816751990457;
    ASSERT_NEAR(numerical_Sig_xy_nl_, analytical_Sig_xy_nl_, 1.e-5);

    double numerical_Sig_zx_nl_ = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zx_nl_ = +5.654491186524116;
    ASSERT_NEAR(numerical_Sig_zx_nl_, analytical_Sig_zx_nl_, 1.e-5);

    double numerical_Sig_zy_nl_ = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz);
    double analytical_Sig_zy_nl_ = +6.596906384278137;
    ASSERT_NEAR(numerical_Sig_zy_nl_, analytical_Sig_zy_nl_, 1.e-5);

    x = -1.2; y = 1.4; z = -1.5;
    numerical_Sig_xx_nl_ = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xx_nl_ = -2.0657867576932993;
    ASSERT_NEAR(numerical_Sig_xx_nl_, analytical_Sig_xx_nl_, 1.e-5);

    numerical_Sig_yy_nl_ = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_yy_nl_ = -3.656171942425667;
    ASSERT_NEAR(numerical_Sig_yy_nl_, analytical_Sig_yy_nl_, 1.e-5);

    numerical_Sig_zz_nl_ = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zz_nl_ = 0.543869725910124;
    ASSERT_NEAR(numerical_Sig_zz_nl_, analytical_Sig_zz_nl_, 1.e-5);

    numerical_Sig_xy_nl_ = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xy_nl_ = 5.13816751990457;
    ASSERT_NEAR(numerical_Sig_xy_nl_, analytical_Sig_xy_nl_, 1.e-5);

    numerical_Sig_zx_nl_ = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zx_nl_ = -5.654491186524116;
    ASSERT_NEAR(numerical_Sig_zx_nl_, analytical_Sig_zx_nl_, 1.e-5);

    numerical_Sig_zy_nl_ = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zy_nl_ = +6.596906384278137;
    ASSERT_NEAR(numerical_Sig_zy_nl_, analytical_Sig_zy_nl_, 1.e-5);


    x = -1.2; y = -1.4; z = -1.5;
    numerical_Sig_xx_nl_ = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xx_nl_ = -2.0657867576932993;
    ASSERT_NEAR(numerical_Sig_xx_nl_, analytical_Sig_xx_nl_, 1.e-5);

    numerical_Sig_yy_nl_ = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_yy_nl_ = -3.656171942425667;
    ASSERT_NEAR(numerical_Sig_yy_nl_, analytical_Sig_yy_nl_, 1.e-5);

    numerical_Sig_zz_nl_ = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zz_nl_ = 0.543869725910124;
    ASSERT_NEAR(numerical_Sig_zz_nl_, analytical_Sig_zz_nl_, 1.e-5);

    numerical_Sig_xy_nl_ = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xy_nl_ = -5.13816751990457;
    ASSERT_NEAR(numerical_Sig_xy_nl_, analytical_Sig_xy_nl_, 1.e-5);

    numerical_Sig_zx_nl_ = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zx_nl_ = -5.654491186524116;
    ASSERT_NEAR(numerical_Sig_zx_nl_, analytical_Sig_zx_nl_, 1.e-5);

    numerical_Sig_zy_nl_ = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zy_nl_ = -6.596906384278137;
    ASSERT_NEAR(numerical_Sig_zy_nl_, analytical_Sig_zy_nl_, 1.e-5);

    x = 1.2; y = -1.4; z = -1.5;
    numerical_Sig_xx_nl_ = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xx_nl_ = -2.0657867576932993;
    ASSERT_NEAR(numerical_Sig_xx_nl_, analytical_Sig_xx_nl_, 1.e-5);

    numerical_Sig_yy_nl_ = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_yy_nl_ = -3.656171942425667;
    ASSERT_NEAR(numerical_Sig_yy_nl_, analytical_Sig_yy_nl_, 1.e-5);

    numerical_Sig_zz_nl_ = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zz_nl_ = 0.543869725910124;
    ASSERT_NEAR(numerical_Sig_zz_nl_, analytical_Sig_zz_nl_, 1.e-5);

    numerical_Sig_xy_nl_ = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_xy_nl_ = 5.13816751990457;
    ASSERT_NEAR(numerical_Sig_xy_nl_, analytical_Sig_xy_nl_, 1.e-5);

    numerical_Sig_zx_nl_ = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zx_nl_ = +5.654491186524116;
    ASSERT_NEAR(numerical_Sig_zx_nl_, analytical_Sig_zx_nl_, 1.e-5);

    numerical_Sig_zy_nl_ = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz);
    analytical_Sig_zy_nl_ = -6.596906384278137;
    ASSERT_NEAR(numerical_Sig_zy_nl_, analytical_Sig_zy_nl_, 1.e-5);
}

// shear load

TEST(radial_analytical_negative_Z, Stress_sl_){
    double x = 1.2, y = 1.4, z = -1.5;
    double a = 1.2;
    double G = 0.56, nu = 0.3;
    double px = 155., py = 230.;

    double numerical_Sig_xx_sl_ = bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_xx_sl_ = 6.376791154478333;
    ASSERT_NEAR(numerical_Sig_xx_sl_, analytical_Sig_xx_sl_, 1.e-5);

    double numerical_Sig_yy_sl_ = bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_yy_sl_ = 9.2234549711472;
    ASSERT_NEAR(numerical_Sig_yy_sl_, analytical_Sig_yy_sl_, 1.e-5);

    double numerical_Sig_zz_sl_ = bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zz_sl_ = 17.600989722758893;
    ASSERT_NEAR(numerical_Sig_zz_sl_, analytical_Sig_zz_sl_, 1.e-5);

    double numerical_Sig_xy_sl_ = bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_xy_sl_ = 9.89471480535422;
    ASSERT_NEAR(numerical_Sig_xy_sl_, analytical_Sig_xy_sl_, 1.e-5);

    double numerical_Sig_zx_sl_ = bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zx_sl_ = -11.60516358682423;
    ASSERT_NEAR(numerical_Sig_zx_sl_, analytical_Sig_zx_sl_, 1.e-5);

    double numerical_Sig_zy_sl_ = bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py);
    double analytical_Sig_zy_sl_ = -12.62186684505416;
    ASSERT_NEAR(numerical_Sig_zy_sl_, analytical_Sig_zy_sl_, 1.e-5);

    x = -1.2; y = 1.4; z = -1.5;

    numerical_Sig_xx_sl_ = bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xx_sl_ = 2.560211770054411;
    ASSERT_NEAR(numerical_Sig_xx_sl_, analytical_Sig_xx_sl_, 1.e-5);

    numerical_Sig_yy_sl_ = bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_yy_sl_ = 1.61623207460115;
    ASSERT_NEAR(numerical_Sig_yy_sl_, analytical_Sig_yy_sl_, 1.e-5);

    numerical_Sig_zz_sl_ = bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zz_sl_ = 4.712075988770087;
    ASSERT_NEAR(numerical_Sig_zz_sl_, analytical_Sig_zz_sl_, 1.e-5);

    numerical_Sig_xy_sl_ = bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xy_sl_ = -2.7809967591087665;
    ASSERT_NEAR(numerical_Sig_xy_sl_, analytical_Sig_xy_sl_, 1.e-5);

    numerical_Sig_zx_sl_ = bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zx_sl_ = 6.77367327453034;
    ASSERT_NEAR(numerical_Sig_zx_sl_, analytical_Sig_zx_sl_, 1.e-5);

    numerical_Sig_zy_sl_ = bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zy_sl_ = -0.2361289602282864;
    ASSERT_NEAR(numerical_Sig_zy_sl_, analytical_Sig_zy_sl_, 1.e-5);

    x = 1.2; y = -1.4; z = -1.5;

    numerical_Sig_xx_sl_ = bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xx_sl_ = -2.560211770054411;
    ASSERT_NEAR(numerical_Sig_xx_sl_, analytical_Sig_xx_sl_, 1.e-5);

    numerical_Sig_yy_sl_ = bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_yy_sl_ = -1.61623207460115;
    ASSERT_NEAR(numerical_Sig_yy_sl_, analytical_Sig_yy_sl_, 1.e-5);

    numerical_Sig_zz_sl_ = bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zz_sl_ = -4.712075988770087;
    ASSERT_NEAR(numerical_Sig_zz_sl_, analytical_Sig_zz_sl_, 1.e-5);

    numerical_Sig_xy_sl_ = bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xy_sl_ = 2.7809967591087665;
    ASSERT_NEAR(numerical_Sig_xy_sl_, analytical_Sig_xy_sl_, 1.e-5);

    numerical_Sig_zx_sl_ = bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zx_sl_ = 6.77367327453034;
    ASSERT_NEAR(numerical_Sig_zx_sl_, analytical_Sig_zx_sl_, 1.e-5);

    numerical_Sig_zy_sl_ = bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zy_sl_ = -0.2361289602282864;
    ASSERT_NEAR(numerical_Sig_zy_sl_, analytical_Sig_zy_sl_, 1.e-5);


    x = -1.2; y = -1.4; z = -1.5;

    numerical_Sig_xx_sl_ = bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xx_sl_ = -6.376791154478333;
    ASSERT_NEAR(numerical_Sig_xx_sl_, analytical_Sig_xx_sl_, 1.e-5);

    numerical_Sig_yy_sl_ = bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_yy_sl_ = -9.2234549711472;
    ASSERT_NEAR(numerical_Sig_yy_sl_, analytical_Sig_yy_sl_, 1.e-5);

    numerical_Sig_zz_sl_ = bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zz_sl_ = -17.600989722758893;
    ASSERT_NEAR(numerical_Sig_zz_sl_, analytical_Sig_zz_sl_, 1.e-5);

    numerical_Sig_xy_sl_ = bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_xy_sl_ = -9.89471480535422;
    ASSERT_NEAR(numerical_Sig_xy_sl_, analytical_Sig_xy_sl_, 1.e-5);

    numerical_Sig_zx_sl_ = bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zx_sl_ = -11.60516358682423;
    ASSERT_NEAR(numerical_Sig_zx_sl_, analytical_Sig_zx_sl_, 1.e-5);

    numerical_Sig_zy_sl_ = bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py);
    analytical_Sig_zy_sl_ = -12.62186684505416;
    ASSERT_NEAR(numerical_Sig_zy_sl_, analytical_Sig_zy_sl_, 1.e-5);



}
