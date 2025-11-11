#include <gtest/gtest.h>

#include <iostream>

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/Array2D.h>

#include "elements/triangle.h"
#include "elements/rectangle.h"

#include "elasticity/fullspace_iso_2d_hydrostatic_eigenstrain/hydrostatic_eigenstrain_polygon.hh"

TEST(TwoDT0S0_V, test1) {

    // test tolerance 
    double tol = 1e-4;

    // Column major !!
    il::Array2D<double> tri_vertices{il::value, {
        {1.0,       -0.5,       -0.5},
        {0.0,       0.8660254,  -0.8660254},
        {0.0,       0.0,        0.0}
    }};

    bigwham::Triangle<0> triangle;
    triangle.SetElement(tri_vertices);

    double nu = 0.25;
    double G = 1.0;

    // ---------- Outside inclusion
    il::StaticArray<double, 2> xy_obs{il::value, {-1.0, 0.0}};
    il::StaticArray<double, 2> n_obs{il::value, {1.0, 0.0}};

    auto t_i = bigwham::V_twoD_polygon_0(
        triangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.41189862, tol) << "Point outside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point outside inclusion : traction mismatch on component 1";

    // ---------- On inclusion boundary
    xy_obs[0] = -0.5;
    xy_obs[1] = 0.;

    t_i = bigwham::V_twoD_polygon_0(
        triangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.79810112, tol) << "Point on inclusion boundary : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point on inclusion boundary : traction mismatch on component 1";


    // ---------- Inside inclusion
    xy_obs[0] = 0.0;
    xy_obs[1] = 0.0;

    t_i = bigwham::V_twoD_polygon_0(
        triangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -1.33333334, tol) << "Point inside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point inside inclusion : traction mismatch on component 1";
}

TEST(TwoDR0S0_V, test1) {

    // test tolerance 
    double tol = 1e-4;

    // Column major !!
    il::Array2D<double> rec_vertices{il::value, {
        {0.0,   1.0,    1.0,    0.0},
        {0.0,   0.0,    1.0,    1.0},
        {0.0,   0.0,    0.0,    0.0}
    }};     
    
    bigwham::Rectangle<0> rectangle;
    rectangle.SetElement(rec_vertices);

    double nu = 0.25;
    double G = 1.0;

    // ---------- Outside inclusion
    il::StaticArray<double, 2> xy_obs{il::value, {1.5, 0.5}};
    il::StaticArray<double, 2> n_obs{il::value, {1.0, 0.0}};

    auto t_i = bigwham::V_twoD_polygon_0(
        rectangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.39355631, tol) << "Point outside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point outside inclusion : traction mismatch on component 1";

    // ---------- On inclusion boundary
    xy_obs[0] = 1.0;
    xy_obs[1] = 0.5;

    t_i = bigwham::V_twoD_polygon_0(
        rectangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.93977702, tol) << "Point on inclusion boundary : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point on inclusion boundary : traction mismatch on component 1";


    // ---------- Inside inclusion
    xy_obs[0] = 0.5;
    xy_obs[1] = 0.5;

    t_i = bigwham::V_twoD_polygon_0(
        rectangle,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -1.33333334, tol) << "Point inside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point inside inclusion : traction mismatch on component 1";
}