#include <gtest/gtest.h>

#include <iostream>

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/Array2D.h>

#include "elements/triangle.h"
#include "elements/rectangle.h"

#include "elasticity/fullspace_iso_3daxis_hydrostatic_eigenstrain/hydrostatic_eigenstrain_3daxis.hh"

TEST(ThreeDAxisR0S0_V, test1) {

    // test tolerance 
    double tol = 1e-4;

    // Column major !!
    il::Array2D<double> rec_vertices{il::value, {
        {1.0,    2.0,  2.0,  1.0},
        {-0.5,  -0.5,  0.5,  0.5},
        {0.0,    0.0,  0.0,  0.0}
    }};     
    
    bigwham::Rectangle<0> rectangle;
    rectangle.SetElement(rec_vertices);

    double nu = 0.33;
    double E = 1.0;
    double G = E / ( 2 * (1 + nu));

    double phi;
    il::StaticArray<double, 2> rz_obs;
    il::StaticArray<double, 2> n_obs;
    std::array<double, 4> e_ij;
    il::StaticArray<double, 2> t_i;

    // ---------- Point (R, Z) = (0.5, 0.0)
    rz_obs = {il::value, {0.5, 0.0}};

    phi = phi_3daxis(rectangle, rz_obs);
    ASSERT_NEAR(phi, -6.3298613193309956, tol) << "Point (0.5, 0.0), phi incorrect";

    e_ij = strain_3daxis(rectangle, rz_obs, nu);
    ASSERT_NEAR(e_ij[0], 2.7613091275549917e-01, tol) << "Point (0.5, 0.0), err incorrect";
    ASSERT_NEAR(e_ij[1], -5.0253876889295490e-01, tol) << "Point (0.5, 0.0), ezz incorrect";
    ASSERT_NEAR(e_ij[2], 2.2640782188759381e-01, tol) << "Point (0.5, 0.0), ett incorrect";
    ASSERT_NEAR(e_ij[3], 0.0, tol) << "Point (0.5, 0.0), erz incorrect";


    // ---------- Point (R, Z) = (1.5, 0.0)
    rz_obs = {il::value, {1.5, 0.0}};

    phi = phi_3daxis(rectangle, rz_obs);
    ASSERT_NEAR(phi, -7.0404287362230731, tol) << "Point (1.5, 0.0), phi incorrect";

    e_ij = strain_3daxis(rectangle, rz_obs, nu);
    ASSERT_NEAR(e_ij[0], -8.2937320144436533e-01, tol) << "Point (1.5, 0.0), err incorrect";
    ASSERT_NEAR(e_ij[1], -1.0146720891433476e+00, tol) << "Point (1.5, 0.0), ezz incorrect";
    ASSERT_NEAR(e_ij[2], -1.4102930389908538e-01, tol) << "Point (1.5, 0.0), ett incorrect";
    ASSERT_NEAR(e_ij[3], 0.0, tol) << "Point (1.5, 0.0), erz incorrect";

    // ---------- Point (R, Z) = (1.5, 1.0)
    rz_obs = {il::value, {1.5, 1.0}};

    phi = phi_3daxis(rectangle, rz_obs);
    ASSERT_NEAR(phi, -4.8577583346682047, tol) << "Point (1.5, 1.0), phi incorrect";

    e_ij = strain_3daxis(rectangle, rz_obs, nu);
    ASSERT_NEAR(e_ij[0], -2.0323619846676805e-01, tol) << "Point (1.5, 1.0), err incorrect";
    ASSERT_NEAR(e_ij[1], 2.9434048443492816e-01, tol) << "Point (1.5, 1.0), ezz incorrect";
    ASSERT_NEAR(e_ij[2], -9.1104201949739824e-02, tol) << "Point (1.5, 1.0), ett incorrect";
    ASSERT_NEAR(e_ij[3], 9.1737853285520185e-02, tol) << "Point (1.5, 1.0), erz incorrect";

    // ---------- Point (R, Z) = (2.0, 0.0), Normal n = [1. 0.]
    rz_obs = {il::value, {2.0, 0.0}};
    n_obs = {il::value, {1.0, 0.0}};

    t_i = V_threeDAxis_polygon_0(rectangle, rz_obs, n_obs, G, nu);
    ASSERT_NEAR(t_i[0], 7.2701628798701257e-01, tol) << "Point (R, Z) = (2.0, 0.0), Normal n = [1. 0.], t_r incorrect";
     ASSERT_NEAR(t_i[1], 0.0, tol) << "Point (R, Z) = (2.0, 0.0), Normal n = [1. 0.], t_z incorrect";

    // ---------- Point (R, Z) = (1.5, 0.5), Normal n = [0. 1.]
    rz_obs = {il::value, {1.5, 0.5}};
    n_obs = {il::value, {0.0, 1.0}};

    t_i = V_threeDAxis_polygon_0(rectangle, rz_obs, n_obs, G, nu);
    ASSERT_NEAR(t_i[0], 6.9342050469404032e-02, tol) << "Point (R, Z) = (1.5, 0.5), Normal n = [0. 1.], t_r incorrect";
    ASSERT_NEAR(t_i[1], 5.1622626201541899e-01, tol) << "Point (R, Z) = (1.5, 0.5), Normal n = [0. 1.], t_z incorrect";
}