#include <gtest/gtest.h>

#include <iostream>

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

#include "elasticity/fullspace_iso_2d_hydrostatic_eigenstrain_triangle_segment/elastic_2dT0S0_V_element.hh"

TEST(TwoDT0S0_V, test1) {

    // test tolerance 
    double tol = 1e-4;

    // Unit triangle 
    std::array<std::array<double, 2>, 3> tri_vertices{{
        {1.0,       0.0},
        {-0.5,      0.8660254},
        {-0.5,     -0.8660254}
    }};

    double nu = 0.25;
    double G = 1.0;

    // ---------- Outside inclusion
    il::StaticArray<double, 2> xy_obs{il::value, {-1.0, 0.0}};
    il::StaticArray<double, 2> n_obs{il::value, {1.0, 0.0}};
    std::cout << "Point outside inclusion : xy =(" << xy_obs[0] << ", " << xy_obs[1] << ")\n";

    auto t_i = bigwham::V_twoD_triangle_0(
        tri_vertices,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.41189862, tol) << "Point outside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point outside inclusion : traction mismatch on component 1";

    // ---------- On inclusion boundary
    xy_obs[0] = -0.5;
    xy_obs[1] = 0.;
    std::cout << "Point on inclusion boundary : xy =(" << xy_obs[0] << ", " << xy_obs[1] << ")\n";

    t_i = bigwham::V_twoD_triangle_0(
        tri_vertices,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -0.79810112, tol) << "Point on inclusion boundary : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point on inclusion boundary : traction mismatch on component 1";


    // ---------- Inside inclusion
    xy_obs[0] = 0.0;
    xy_obs[1] = 0.0;
    std::cout << "Point inside inclusion : xy =(" << xy_obs[0] << ", " << xy_obs[1] << ")\n";

    t_i = bigwham::V_twoD_triangle_0(
        tri_vertices,
        xy_obs,
        n_obs,
        G, nu
    );

    ASSERT_NEAR(t_i[0], -1.33333334, tol) << "Point inside inclusion : traction mismatch on component 0";
    ASSERT_NEAR(t_i[1], 0., tol) << "Point inside inclusion : traction mismatch on component 1";
}