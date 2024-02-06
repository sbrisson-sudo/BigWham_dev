//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>
#include <gtest/gtest.h>
#include "elastic3DR0_element_benchmark.h"
#include "elasticity/_test/penny_shaped_crack_analytical_sol.h"
#include <BigWhamIO.h>

//--------------------------------------------------------------------------


TEST(R0, benchmark1_displ) {
    /*
     * This benchmark tests the kernel 3DR0.
     * Consider a rectangular crack lying on a plane xy
     * The edges of such crack are 2a (a=2 )and 2b (b=4) respectively parallel to the axes x and y.
     *
     *  We impose a local displacement discontinuity of {1,1,1}
     *  We then check the stress tensor at {1., 2., 3.} and at {-1., -2., -3.}
     *  The crack plane is discretized with 8x8 rectangular elements.
     */
    const int max_leaf_size = 100; //max size of the smallest sub block
    const double eta = 0.;// it governs the severity of the sub blocks partitioning based on the distance
    const double eps_aca = 0.0001; // governs the stopping criterion for the approximation
    Rect_mesh_1 mymesh; // mesh of size: 4 in dir x
    //               8 in dir y
    // mesh centered at (0,0,0) with normal (0,0,1)

    const std::vector<double> properties = {100, 0.2}; // Young Modulus , Poisson's ratio

    /* --------------------------------------------------
     * create traction HMAT to test the stresses
     * Benchmark 1
     * --------------------------------------------------
     */
    const std::string tractionKernel = "3DR0";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points
    std::vector<double> dd_solution(mymesh.nelts * 3), obsPoints = {1., 2., 3., -1., -2., -3.};

    for (int x = 0; x < mymesh.nelts * 3; ++x) dd_solution[x] = 1.;

    std::vector<double> stressATpoints = tractionHMAT.computeStresses(dd_solution, obsPoints, 2, properties,
                                                                      mymesh.coor, mymesh.conn, true);

    // solution point {1,2,3}:
    double x, y, z, a, b, G, nu;
    x = 1. - 0.;
    y = 2. - 0.;
    z = 3. - 0.;
    a = 2.;
    b = 4.;
    G = properties[0] / (2.0 * (1 + properties[1]));
    nu = properties[1];

    il::StaticArray2D<double, 3, 6> Stress;
    Stress = bie::StressesKernelR0(x, y, z, a, b, G, nu);

    std::vector<double> stress_analytical(6, 0.);
    stress_analytical[0] = Stress(0, 0) + Stress(1, 0) + Stress(2, 0); // sxx
    stress_analytical[1] = Stress(0, 1) + Stress(1, 1) + Stress(2, 1); // syy
    stress_analytical[2] = Stress(0, 2) + Stress(1, 2) + Stress(2, 2); // szz
    stress_analytical[3] = Stress(0, 3) + Stress(1, 3) + Stress(2, 3); // sxy
    stress_analytical[4] = Stress(0, 4) + Stress(1, 4) + Stress(2, 4); // sxz
    stress_analytical[5] = Stress(0, 5) + Stress(1, 5) + Stress(2, 5); // syz

    ASSERT_NEAR(stress_analytical[0], stressATpoints[0], 1.e-5);
    ASSERT_NEAR(stress_analytical[1], stressATpoints[1], 1.e-5);
    ASSERT_NEAR(stress_analytical[2], stressATpoints[2], 1.e-5);
    ASSERT_NEAR(stress_analytical[3], stressATpoints[3], 1.e-5);
    ASSERT_NEAR(stress_analytical[4], stressATpoints[4], 1.e-5);
    ASSERT_NEAR(stress_analytical[5], stressATpoints[5], 1.e-5);

    // solution at point {-1.,-2.,-3.}:
    x = -1. - 0.;
    y = -2. - 0.;
    z = -3. - 0.;
    Stress = bie::StressesKernelR0(x, y, z, a, b, G, nu);
    stress_analytical[0] = +Stress(0, 0) + Stress(1, 0) + Stress(2, 0); // sxx
    stress_analytical[1] = +Stress(0, 1) + Stress(1, 1) + Stress(2, 1); // syy
    stress_analytical[2] = +Stress(0, 2) + Stress(1, 2) + Stress(2, 2); // szz
    stress_analytical[3] = +Stress(0, 3) + Stress(1, 3) + Stress(2, 3); // sxy
    stress_analytical[4] = +Stress(0, 4) + Stress(1, 4) + Stress(2, 4); // sxz
    stress_analytical[5] = +Stress(0, 5) + Stress(1, 5) + Stress(2, 5); // syz

    ASSERT_NEAR(stress_analytical[0], stressATpoints[6], 1.e-5);
    ASSERT_NEAR(stress_analytical[1], stressATpoints[7], 1.e-5);
    ASSERT_NEAR(stress_analytical[2], stressATpoints[8], 1.e-5);
    ASSERT_NEAR(stress_analytical[3], stressATpoints[9], 1.e-5);
    ASSERT_NEAR(stress_analytical[4], stressATpoints[10], 1.e-5);
    ASSERT_NEAR(stress_analytical[5], stressATpoints[11], 1.e-5);

    /* --------------------------------------------------
     * create displacement HMAT to test the displacements
     * Benchmark 1
     * --------------------------------------------------
     */
    const std::string displacementKernel = "3DR0_displ";
    Bigwhamio displacementHMAT;
    displacementHMAT.set(mymesh.coor,mymesh.conn, displacementKernel, properties, max_leaf_size, eta, eps_aca);

    // get displacements at points {1., 2., 3., -1., -2., -3.}

    for (int i = 0; i < mymesh.nelts * 3; ++i) dd_solution[i] = 1.;

    std::vector<double> displATpoints = displacementHMAT.computeDisplacements(dd_solution, obsPoints, 2, properties,
                                                                                  mymesh.coor, mymesh.conn, true);

    // solution at point {1.,2.,3.}:
    x = 1. - 0.;
    y = 2. - 0.;
    z = 3. - 0.;
    il::Array2D<double> Displ{3,3};
    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
    //   0      -> |       Ux,            Ux,             Ux            |
    //   1      -> |       Uy,            Uy,             Uy            |
    //   2      -> |       Uz,            Uz,             Uz            |

    std::vector<double> displ_analytical(3, 0.);
    displ_analytical[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    ASSERT_NEAR(displ_analytical[0], displATpoints[0], 1.e-5);
    ASSERT_NEAR(displ_analytical[1], displATpoints[1], 1.e-5);
    ASSERT_NEAR(displ_analytical[2], displATpoints[2], 1.e-5);

    // solution at point {-1.,-2.,-3.}:
    x = -1. - 0.;
    y = -2. - 0.;
    z = -3. - 0.;
    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    displ_analytical[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    ASSERT_NEAR(displ_analytical[0], displATpoints[3], 1.e-5);
    ASSERT_NEAR(displ_analytical[1], displATpoints[4], 1.e-5);
    ASSERT_NEAR(displ_analytical[2], displATpoints[5], 1.e-5);
}

TEST(R0, benchmark2_displ) {
    /*
     * This benchmark tests the kernel 3DR0.
     * Consider a rectangular crack initially lying on a plane xy
     * The edges of such crack are 2a (a=2 )and 2b (b=4) respectively initially parallel to the axes x and y.
     *
     *  Finally apply a translation of a vector {10.,20.,30.}
     *
     *  We impose a local displacement discontinuity of {1,1,1}
     *  We then check the stress tensor at {1., 2., 3.} and at {-1., -2., -3.}
     *  The crack plane is discretized with 8x8 rectangular elements.
     */
    const int max_leaf_size = 1000; //max size of the smallest sub block
    const double eta = 0.;// it governs the severity of the sub blocks partitioning based on the distance
    const double eps_aca = 0.0001; // governs the stopping criterion for the approximation
    Rect_mesh_2 mymesh; // mesh of size: 4 in dir x
    //               8 in dir y
    // mesh centered at (0,0,0) with normal (0,0,1)

    const std::vector<double> properties = {100, 0.2}; // Young Modulus , Poisson's ratio

    /* --------------------------------------------------
     * create traction HMAT to test the stresses
     * Benchmark 2
     * --------------------------------------------------
     */
    const std::string tractionKernel = "3DR0";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points
    std::vector<double> dd_solution(mymesh.nelts * 3), obsPoints = {1., 2., 3., -1.,-2.,-3.};

    for (int x = 0; x < mymesh.nelts * 3; ++x) dd_solution[x] = 1.;

    std::vector<double> stressATpoints = tractionHMAT.computeStresses(dd_solution, obsPoints, 2, properties,
                                                                      mymesh.coor, mymesh.conn, true);

    // solution at point : {1., 2., 3.}
    double x, y, z, a, b, G, nu;
    x = +1. - 10.;
    y = +2. - 20.;
    z = +3. - 30.;
    a = 2.;
    b = 4.;
    G = properties[0] / (2.0 * (1 + properties[1]));
    nu = properties[1];

    il::StaticArray2D<double, 3, 6> Stress;
    Stress = bie::StressesKernelR0(x, y, z, a, b, G, nu);

    std::vector<double> stress_analytical(6, 0.);
    stress_analytical[0] = Stress(0, 0) + Stress(1, 0) + Stress(2, 0); // sxx
    stress_analytical[1] = Stress(0, 1) + Stress(1, 1) + Stress(2, 1); // syy
    stress_analytical[2] = Stress(0, 2) + Stress(1, 2) + Stress(2, 2); // szz
    stress_analytical[3] = Stress(0, 3) + Stress(1, 3) + Stress(2, 3); // sxy
    stress_analytical[4] = Stress(0, 4) + Stress(1, 4) + Stress(2, 4); // sxz
    stress_analytical[5] = Stress(0, 5) + Stress(1, 5) + Stress(2, 5); // syz

    ASSERT_NEAR(stress_analytical[0], stressATpoints[0], 1.e-5);
    ASSERT_NEAR(stress_analytical[1], stressATpoints[1], 1.e-5);
    ASSERT_NEAR(stress_analytical[2], stressATpoints[2], 1.e-5);
    ASSERT_NEAR(stress_analytical[3], stressATpoints[3], 1.e-5);
    ASSERT_NEAR(stress_analytical[4], stressATpoints[4], 1.e-5);
    ASSERT_NEAR(stress_analytical[5], stressATpoints[5], 1.e-5);

    // solution at point : {-1., -2., -3.}
    x = -1. - 10.;
    y = -2. - 20.;
    z = -3. - 30.;

    Stress = bie::StressesKernelR0(x, y, z, a, b, G, nu);

    stress_analytical[0] = Stress(0, 0) + Stress(1, 0) + Stress(2, 0); // sxx
    stress_analytical[1] = Stress(0, 1) + Stress(1, 1) + Stress(2, 1); // syy
    stress_analytical[2] = Stress(0, 2) + Stress(1, 2) + Stress(2, 2); // szz
    stress_analytical[3] = Stress(0, 3) + Stress(1, 3) + Stress(2, 3); // sxy
    stress_analytical[4] = Stress(0, 4) + Stress(1, 4) + Stress(2, 4); // sxz
    stress_analytical[5] = Stress(0, 5) + Stress(1, 5) + Stress(2, 5); // syz

    ASSERT_NEAR(stress_analytical[0], stressATpoints[6], 1.e-5);
    ASSERT_NEAR(stress_analytical[1], stressATpoints[7], 1.e-5);
    ASSERT_NEAR(stress_analytical[2], stressATpoints[8], 1.e-5);
    ASSERT_NEAR(stress_analytical[3], stressATpoints[9], 1.e-5);
    ASSERT_NEAR(stress_analytical[4], stressATpoints[10], 1.e-5);
    ASSERT_NEAR(stress_analytical[5], stressATpoints[11], 1.e-5);

    /* --------------------------------------------------
     * create displacement HMAT to test the displacements
     * Benchmark 2
     * --------------------------------------------------
     */

    const std::string displacementKernel = "3DR0_displ";
    Bigwhamio displacementHMAT;
    displacementHMAT.set(mymesh.coor,mymesh.conn, displacementKernel, properties, max_leaf_size, eta, eps_aca);

    // get displacements at points {1., 2., 3., -1., -2., -3.}

    for (int i = 0; i < mymesh.nelts * 3; ++i) dd_solution[i] = 1.;

    std::vector<double> displATpoints = displacementHMAT.computeDisplacements(dd_solution, obsPoints, 2, properties,
                                                                              mymesh.coor, mymesh.conn, true);

    // solution at point {1.,2.,3.}:
    x = +1. - 10.;
    y = +2. - 20.;
    z = +3. - 30.;
    il::Array2D<double> Displ{3,3};
    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
    //   0      -> |       Ux,            Ux,             Ux            |
    //   1      -> |       Uy,            Uy,             Uy            |
    //   2      -> |       Uz,            Uz,             Uz            |

    std::vector<double> displ_analytical(3, 0.);
    displ_analytical[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    ASSERT_NEAR(displ_analytical[0], displATpoints[0], 1.e-5);
    ASSERT_NEAR(displ_analytical[1], displATpoints[1], 1.e-5);
    ASSERT_NEAR(displ_analytical[2], displATpoints[2], 1.e-5);

    // solution at point {-1.,-2.,-3.}:
    x = -1. - 10.;
    y = -2. - 20.;
    z = -3. - 30.;
    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    displ_analytical[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    ASSERT_NEAR(displ_analytical[0], displATpoints[3], 1.e-5);
    ASSERT_NEAR(displ_analytical[1], displATpoints[4], 1.e-5);
    ASSERT_NEAR(displ_analytical[2], displATpoints[5], 1.e-5);
}

TEST(R0, benchmark3_displ) {
    /*
     * This benchmark tests the kernel 3DR0.
     * Consider a rectangular crack initially lying on a plane xy
     * The edges of such crack are 2a (a=2 )and 2b (b=4) respectively initially parallel to the axes x and y.
     * Then we apply a rotation to te plane according to the following matrix: (local to global matrix)
     *
     *  |-0.690578, 0.597022, -0.408248|
     *  |0.568711, 0.0995037, -0.816497|
     *  |-0.446844, -0.79603, -0.408248|
     *
     *  Finally apply a translation of a vector {10.,20.,30.}
     *
     *  We impose a local displacement discontinuity of {1,1,1}
     *  We then check the stress tensor at the origin {0.,0.,0.} and {100., 100., 100.}
     *  The crack plane is discretized with 8x8 rectangular elements.
     */
    const int max_leaf_size = 1000; //max size of the smallest sub block
    const double eta = 0.;// it governs the severity of the sub blocks partitioning based on the distance
    const double eps_aca = 0.00001; // governs the stopping criterion for the approximation
    Rect_mesh_3 mymesh; // mesh of size: 4 in dir x
    //               8 in dir y
    // mesh centered at (0,0,0) with normal (0,0,1)

    const std::vector<double> properties = {100, 0.2}; // Young Modulus , Poisson's ratio

    /* --------------------------------------------------
     * create traction HMAT to test the stresses
     * Benchmark 3
     * --------------------------------------------------
     */
    const std::string tractionKernel = "3DR0";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points by prescribing global DD
    std::vector<double>  dd_solution_g(mymesh.nelts * 3), obsPoints = {0.,0.,0., 100., 100., 100.};

    for(int x = 0; x < mymesh.nelts * 3; x=x+3)
    { dd_solution_g[x] = -0.5018039999999999;
      dd_solution_g[x+1] = -0.1482823;
      dd_solution_g[x+2] = -1.651122;}

    std::vector<double> stressATpoints_0 = tractionHMAT.computeStresses(dd_solution_g, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, true);

    // get stress at points by prescribing local DD
    std::vector<double>  dd_solution_l(mymesh.nelts * 3);

    for(int x = 0; x < mymesh.nelts * 3; x=x+3)
    { dd_solution_l[x] = 1;
      dd_solution_l[x+1] = 1;
      dd_solution_l[x+2] = 1;}

    std::vector<double> stressATpoints = tractionHMAT.computeStresses(dd_solution_l, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, false);

    // check that, if properly done, prescribing global DD or local DDs does not matter
    ASSERT_NEAR(stressATpoints_0[0], stressATpoints[0], 4.e-4);
    ASSERT_NEAR(stressATpoints_0[1], stressATpoints[1], 4.e-4);
    ASSERT_NEAR(stressATpoints_0[2], stressATpoints[2], 4.e-4);
    ASSERT_NEAR(stressATpoints_0[3], stressATpoints[3], 4.e-4);
    ASSERT_NEAR(stressATpoints_0[4], stressATpoints[4], 4.e-4);
    ASSERT_NEAR(stressATpoints_0[5], stressATpoints[5], 4.e-4);


    // solution at point {0,0,0}:

    // Rotation matricies:
    il::Array2D<double> RT{3,3,0.}; // global to local ...R(g->l)
    RT(0,0) = -0.690578; RT(0,1) = 0.597022; RT(0,2) = -0.408248;
    RT(1,0) = 0.568711;  RT(1,1) = 0.0995037; RT(1,2) = -0.816497;
    RT(2,0) = -0.446844; RT(2,1) = -0.79603; RT(2,2) = -0.408248;

    il::Array2D<double> R{3,3,0.}; // local to global ...R(l->g)
    R(0,0) = -0.690578; R(0,1) = 0.568711;  R(0,2) = -0.446844;
    R(1,0) =  0.597022; R(1,1) = 0.0995037; R(1,2) = -0.79603;
    R(2,0) =  -0.408248; R(2,1) = -0.816497;  R(2,2) =  -0.408248;

    double x, y, z, a, b, G, nu;
    il::Array<double> center_translation ={il::value,{10,20,30}};
    il::Array<double> relative_distance{3,0.};
    relative_distance[0] = 0. - center_translation[0];
    relative_distance[1] = 0. - center_translation[1];
    relative_distance[2] = 0. - center_translation[2];
    relative_distance = il::dot(R,relative_distance);
    x = relative_distance[0]; y = relative_distance[1]; z = relative_distance[2];

    a = 2.; b = 4.;
    G = properties[0] / (2.0 * (1 + properties[1])); nu = properties[1];

    il::StaticArray2D<double, 3, 6> Stress;
    Stress = bie::StressesKernelR0(x, y,  z,  a,  b, G, nu);

    std::vector<double> stress_analytical_l_l(6,0.), stress_analytical_g_g(6,0.);
    stress_analytical_l_l[0] = Stress(0,0) + Stress(1,0) + Stress(2,0); // sxx
    stress_analytical_l_l[1] = Stress(0,1) + Stress(1,1) + Stress(2,1); // syy
    stress_analytical_l_l[2] = Stress(0,2) + Stress(1,2) + Stress(2,2); // szz
    stress_analytical_l_l[3] = Stress(0,3) + Stress(1,3) + Stress(2,3); // sxy
    stress_analytical_l_l[4] = Stress(0,4) + Stress(1,4) + Stress(2,4); // sxz
    stress_analytical_l_l[5] = Stress(0,5) + Stress(1,5) + Stress(2,5); // syz

    il::Array2D<double> StressTensor_l_l{3,3};
    StressTensor_l_l(0,0) = stress_analytical_l_l[0]; // sxx
    StressTensor_l_l(0,1) = stress_analytical_l_l[3]; // sxy
    StressTensor_l_l(0,2) = stress_analytical_l_l[4]; // sxz
    StressTensor_l_l(1,0) = stress_analytical_l_l[3]; // sxy
    StressTensor_l_l(1,1) = stress_analytical_l_l[1]; // syy
    StressTensor_l_l(1,2) = stress_analytical_l_l[5]; // syz
    StressTensor_l_l(2,0) = stress_analytical_l_l[4]; // sxz
    StressTensor_l_l(2,1) = stress_analytical_l_l[5]; // syz
    StressTensor_l_l(2,2) = stress_analytical_l_l[2]; // szz

    //
    //    g = global ref. system
    //    l = local ref. system
    //    Sig = stress tensor
    //
    //    R(l->g)*Sig(l->l)*R(g->l) n(g) = R(l->g) n(l) = n(g)
    //
    il::Array2D<double> StressTensor_g_g =  il::dot(il::dot(RT, StressTensor_l_l ),R );

    stress_analytical_g_g[0] = StressTensor_g_g(0,0) ; // sxx
    stress_analytical_g_g[1] = StressTensor_g_g(1,1) ; // syy
    stress_analytical_g_g[2] = StressTensor_g_g(2,2) ;  // szz
    stress_analytical_g_g[3] = StressTensor_g_g(0,1) ; // sxy
    stress_analytical_g_g[4] = StressTensor_g_g(0,2) ; // sxz
    stress_analytical_g_g[5] = StressTensor_g_g(1,2) ; // syz

    ASSERT_NEAR(stress_analytical_g_g[0], stressATpoints[0], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[1], stressATpoints[1], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[2], stressATpoints[2], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[3], stressATpoints[3], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[4], stressATpoints[4], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[5], stressATpoints[5], 2.e-5);

    // solution at point {100., 100., 100.}:
    relative_distance[0] = 100. - center_translation[0];
    relative_distance[1] = 100. - center_translation[1];
    relative_distance[2] = 100. - center_translation[2];

    relative_distance = il::dot(R,relative_distance);
    x = relative_distance[0]; y = relative_distance[1]; z = relative_distance[2];

    Stress = bie::StressesKernelR0(x, y,  z,  a,  b, G, nu);

    stress_analytical_l_l[0] = Stress(0,0) + Stress(1,0) + Stress(2,0); // sxx
    stress_analytical_l_l[1] = Stress(0,1) + Stress(1,1) + Stress(2,1); // syy
    stress_analytical_l_l[2] = Stress(0,2) + Stress(1,2) + Stress(2,2); // szz
    stress_analytical_l_l[3] = Stress(0,3) + Stress(1,3) + Stress(2,3); // sxy
    stress_analytical_l_l[4] = Stress(0,4) + Stress(1,4) + Stress(2,4); // sxz
    stress_analytical_l_l[5] = Stress(0,5) + Stress(1,5) + Stress(2,5); // syz

    StressTensor_l_l(0,0) = stress_analytical_l_l[0]; // sxx
    StressTensor_l_l(0,1) = stress_analytical_l_l[3]; // sxy
    StressTensor_l_l(0,2) = stress_analytical_l_l[4]; // sxz
    StressTensor_l_l(1,0) = stress_analytical_l_l[3]; // sxy
    StressTensor_l_l(1,1) = stress_analytical_l_l[1]; // syy
    StressTensor_l_l(1,2) = stress_analytical_l_l[5]; // syz
    StressTensor_l_l(2,0) = stress_analytical_l_l[4]; // sxz
    StressTensor_l_l(2,1) = stress_analytical_l_l[5]; // syz
    StressTensor_l_l(2,2) = stress_analytical_l_l[2]; // szz

    //
    //    g = global ref. system
    //    l = local ref. system
    //    Sig = stress tensor
    //
    //    R(l->g)*Sig(l->l)*R(g->l) n(g) = R(l->g) n(l) = n(g)
    //
    StressTensor_g_g =  il::dot(il::dot(RT, StressTensor_l_l ),R );;

    stress_analytical_g_g[0] = StressTensor_g_g(0,0) ; // sxx
    stress_analytical_g_g[1] = StressTensor_g_g(1,1) ; // syy
    stress_analytical_g_g[2] = StressTensor_g_g(2,2) ;  // szz
    stress_analytical_g_g[3] = StressTensor_g_g(0,1) ; // sxy
    stress_analytical_g_g[4] = StressTensor_g_g(0,2) ; // sxz
    stress_analytical_g_g[5] = StressTensor_g_g(1,2) ; // syz

    ASSERT_NEAR(stress_analytical_g_g[0], stressATpoints[6], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[1], stressATpoints[7], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[2], stressATpoints[8], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[3], stressATpoints[9], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[4], stressATpoints[10], 2.e-5);
    ASSERT_NEAR(stress_analytical_g_g[5], stressATpoints[11], 2.e-5);


    /* --------------------------------------------------
     * create displacement HMAT to test the displacements
     * Benchmark 3
     * --------------------------------------------------
     */

    const std::string displacementKernel = "3DR0_displ";
    Bigwhamio displacementHMAT;
    displacementHMAT.set(mymesh.coor,mymesh.conn, displacementKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points by prescribing global DD
    for(int i = 0; i < mymesh.nelts * 3; i=i+3)
    {   dd_solution_g[i] = -0.5687111245916712;
        dd_solution_g[i+1] = -0.09950371902099886;
        dd_solution_g[i+2] = -1.6329931618554523;}

    std::vector<double> displATpoints_0 = displacementHMAT.computeDisplacements(dd_solution_g, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, true);

    // get stress at points by prescribing local DD
    for(int i = 0; i < mymesh.nelts * 3; i=i+3)
    { dd_solution_l[i] = 1;
        dd_solution_l[i+1] = 1;
        dd_solution_l[i+2] = 1;}

    std::vector<double> displATpoints = displacementHMAT.computeDisplacements(dd_solution_l, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, false);

    // check that, if properly done, prescribing global DD or local DDs does not matter
    ASSERT_NEAR(displATpoints_0[0], displATpoints[0], 2.e-5);
    ASSERT_NEAR(displATpoints_0[1], displATpoints[1], 8.e-5);
    ASSERT_NEAR(displATpoints_0[2], displATpoints[2], 6.e-5);
    ASSERT_NEAR(displATpoints_0[3], displATpoints[3], 2.e-5);
    ASSERT_NEAR(displATpoints_0[4], displATpoints[4], 2.e-5);
    ASSERT_NEAR(displATpoints_0[5], displATpoints[5], 2.e-5);

    // solution at point {0,0,0}:

    center_translation = {il::value,{10.,20.,30.}};
    relative_distance[0] = 0. - center_translation[0];
    relative_distance[1] = 0. - center_translation[1];
    relative_distance[2] = 0. - center_translation[2];
    relative_distance = il::dot(R, relative_distance);
    x = relative_distance[0]; y = relative_distance[1]; z = relative_distance[2];


    il::Array2D<double> Displ{3,3};
    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    // expressed in the reference system of the DD element
    // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
    //   0      -> |       Ux,            Ux,             Ux            |
    //   1      -> |       Uy,            Uy,             Uy            |
    //   2      -> |       Uz,            Uz,             Uz            |

    il::Array<double> displ_analytical_l(3, 0.);
    displ_analytical_l[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical_l[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical_l[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    il::Array<double> displ_analytical_g =  il::dot(RT, displ_analytical_l );

    ASSERT_NEAR(displ_analytical_g[0], displATpoints[0], 1.e-5);
    ASSERT_NEAR(displ_analytical_g[1], displATpoints[1], 1.e-5);
    ASSERT_NEAR(displ_analytical_g[2], displATpoints[2], 1.e-5);

    // solution at point {100., 100., 100.}:
    relative_distance[0] = 100. - center_translation[0];
    relative_distance[1] = 100. - center_translation[1];
    relative_distance[2] = 100. - center_translation[2];

    relative_distance = il::dot(R,relative_distance);
    x = relative_distance[0]; y = relative_distance[1]; z = relative_distance[2];

    Displ = bie::DisplacementKernelR0(x, y, z, a, b, nu);
    displ_analytical_l[0] = Displ(0, 0) + Displ(0, 1) + Displ(0, 2); // ux
    displ_analytical_l[1] = Displ(1, 0) + Displ(1, 1) + Displ(1, 2); // uy
    displ_analytical_l[2] = Displ(2, 0) + Displ(2, 1) + Displ(2, 2); // uz

    displ_analytical_g =  il::dot(RT, displ_analytical_l );

    ASSERT_NEAR(displ_analytical_g[0], displATpoints[3], 1.e-5);
    ASSERT_NEAR(displ_analytical_g[1], displATpoints[4], 1.e-5);
    ASSERT_NEAR(displ_analytical_g[2], displATpoints[5], 1.e-5);

}

// the test below is desactivated as we now disabled C++ GMRES
//TEST(R0, benchmark_penny_shaped_crack) {
//    /*
//     * This benchmark tests the kernel 3DR0.
//     * Consider a circular crack initially lying on a plane xy
//     * The crack has radius a=1.2
//     *
//     *  1) We impose an uniform and local load of {0., 0., 160.} over the faces
//     *  2) We impose an uniform and local load of {155., 230., 0.} over the faces
//     *  We then check the stress tensor and the displacements at:
//     *  {1.2, 1.4, +1.5}
//     *  {1.2, 1.4, -1.5}
//     *  The crack plane is discretized with 533 rectangular elements in total.
//     */
//
//    const int max_leaf_size = 1000; //max size of the smallest sub block
//    const double eta = 0.;// it governs the severity of the sub blocks partitioning based on the distance
//    const double eps_aca = 0.00001; // governs the stopping criterion for the approximation
//
//    // mesh
//    Circ_mesh_6 mymesh; // mesh centered at (0,0,0) with normal (0,0,1)
//
//    // elastic properties
//    double G = 0.56;
//    double nu = 0.3;
//    const std::vector<double> properties = {2 * G * (1 + nu), nu}; // Young Modulus , Poisson's ratio
//
//    // kernel
//    const std::string tractionKernel = "3DR0";
//
//    // building the kernel
//    Bigwhamio tractionHMAT;
//    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);
//
//    // get the mesh object to set the rhs
//    il::Array2D<il::int_t> Conn = tractionHMAT.getConn(mymesh.conn);
//    il::Array2D<double> Coor = tractionHMAT.getCoor(mymesh.coor);
//    int p = 0 ; // interpolation order is 0 for the rectangular element
//    bie::Mesh3D mesh3d(Coor, Conn, 0, false);
//
//    // rhs
//    auto * ff = (new  double [mymesh.nelts * 3]); // remember: 3 unknowns per element
//
//    // lhs
//    auto * x0 = (new  double [mymesh.nelts * 3]); // remember: 3 unknowns per element
//
//    // load
//    double px = 0., py = 0., pz = 160.;
//    il::Array<double> load_g{il::value, {-px,-py,-pz}};
//    il::Array<double> load_l{3,0.};
//    int elem = 0;
//    for(int i = 0; i < mymesh.nelts * 3; i= i + 3)
//    {
//        bie::FaceData face = mesh3d.getElementData(elem);
//        il::Array2D<double> R = face.rotationMatrix(true); // g->l rotation matrix
//        load_l = il::dot(R,load_g); // get the local load
//
//        ff[i] =   load_l[0];
//        ff[i+1] = load_l[1];
//        ff[i+2] = load_l[2];
//
//
//        x0[i]=0.;
//        x0[i+1]=0.;
//        x0[i+2]=0.;
//
//        // incremement the element iterator
//        elem ++;
//    }
//
//    // gmres via bigwham object.
//    int nIts = 10000; double tolerance=1.e-6;
//    int stat= tractionHMAT.h_IterativeSolve(ff, false, 1000, il::io, nIts, tolerance, x0);
//
//    std::cout << " Gmress success (0 ?) ";
//    if(stat == 0){std::cout << " YES "; }
//    else { std::cout << " NO ";}
//    std::cout <<"\n";
//
//    std::vector<double> x0_(mymesh.nelts * 3);
//    for (int i=1;i<mymesh.nelts * 3;i++ ){
//        x0_[i]=x0[i];
//    }
//
//    // get stress at obs points
//    std::vector<double>  obsPoints = {1.6, 1.4, +1.5,
//                                      1.6, 1.4, -1.5};
//
//    std::vector<double> stressATpoints = tractionHMAT.computeStresses(x0_, obsPoints, 2, properties , mymesh.coor, mymesh.conn, false);
//
//    double a = 1.2; // crack radius
//    double x = obsPoints[0], y = obsPoints[1], z = obsPoints[2];
//
//    std::vector<double> stress_analytical_(12);
//    stress_analytical_[0] = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    stress_analytical_[1] = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    stress_analytical_[2] = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//    stress_analytical_[3] = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py) ; // sxy
//    stress_analytical_[4] = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py) ; // sxz
//    stress_analytical_[5] = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py) ; // syz
//
//    ASSERT_NEAR(stress_analytical_[0], stressATpoints[0], abs(0.1 * stress_analytical_[0]));
//    ASSERT_NEAR(stress_analytical_[1], stressATpoints[1], abs(0.1 * stress_analytical_[1]));
//    ASSERT_NEAR(stress_analytical_[2], stressATpoints[2], abs(0.1 * stress_analytical_[2]));
//    ASSERT_NEAR(stress_analytical_[3], stressATpoints[3], abs(0.1 * stress_analytical_[3]));
//    ASSERT_NEAR(stress_analytical_[4], stressATpoints[4], abs(0.1 * stress_analytical_[4]));
//    ASSERT_NEAR(stress_analytical_[5], stressATpoints[5], abs(0.1 * stress_analytical_[5]));
//
//    x = obsPoints[3], y = obsPoints[4], z = obsPoints[5];
//    stress_analytical_[6]  = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    stress_analytical_[7]  = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    stress_analytical_[8]  = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//    stress_analytical_[9]  = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py) ; // sxy
//    stress_analytical_[10] = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py) ; // sxz
//    stress_analytical_[11] = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py) ; // syz
//
//    ASSERT_NEAR(stress_analytical_[6], stressATpoints[6], abs(0.1 * stress_analytical_[6]));
//    ASSERT_NEAR(stress_analytical_[7], stressATpoints[7], abs(0.1 * stress_analytical_[7]));
//    ASSERT_NEAR(stress_analytical_[8], stressATpoints[8], abs(0.1 * stress_analytical_[8]));
//    ASSERT_NEAR(stress_analytical_[9], stressATpoints[9], abs(0.1 * stress_analytical_[9]));
//    ASSERT_NEAR(stress_analytical_[10], stressATpoints[10], abs(0.1 * stress_analytical_[10]));
//    ASSERT_NEAR(stress_analytical_[11], stressATpoints[11], abs(0.1 * stress_analytical_[11]));
//
////    std::cout << " \n ********************************* \n";
////    std::cout << " \n relative difference stresses [%]: \n";
////    std::vector<double> relative_err_stress(12);
////    for (int i=0; i<12; i++){
////        relative_err_stress[i]=abs(stress_analytical_[i]-stressATpoints[i])/abs(stress_analytical_[i]);
////        std::cout << relative_err_stress[i] * 100 << "   ";
////    }
//
//    // check displacements
//    std::vector<double> displATpoints = tractionHMAT.computeDisplacements(x0_, obsPoints, 2, properties , mymesh.coor, mymesh.conn, false);
//
//    x = obsPoints[0], y = obsPoints[1], z = obsPoints[2];
//    std::vector<double> displ_analytical_(6);
//    displ_analytical_[0] = bie::Ux_nl_(x, y, z, a, G, nu, pz) + bie::Ux_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    displ_analytical_[1] = bie::Uy_nl_(x, y, z, a, G, nu, pz) + bie::Uy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    displ_analytical_[2] = bie::Uz_nl_(x, y, z, a, G, nu, pz) + bie::Uz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//
//    ASSERT_NEAR(displ_analytical_[0], displATpoints[0], abs(0.1 * displ_analytical_[0]));
//    ASSERT_NEAR(displ_analytical_[1], displATpoints[1], abs(0.1 * displ_analytical_[1]));
//    ASSERT_NEAR(displ_analytical_[2], displATpoints[2], abs(0.1 * displ_analytical_[2]));
//
//    x = obsPoints[3], y = obsPoints[4], z = obsPoints[5];
//    displ_analytical_[3] = bie::Ux_nl_(x, y, z, a, G, nu, pz) + bie::Ux_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    displ_analytical_[4] = bie::Uy_nl_(x, y, z, a, G, nu, pz) + bie::Uy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    displ_analytical_[5] = bie::Uz_nl_(x, y, z, a, G, nu, pz) + bie::Uz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//
//    ASSERT_NEAR(displ_analytical_[3], displATpoints[3], abs(0.1 * displ_analytical_[3]));
//    ASSERT_NEAR(displ_analytical_[4], displATpoints[4], abs(0.1 * displ_analytical_[4]));
//    ASSERT_NEAR(displ_analytical_[5], displATpoints[5], abs(0.1 * displ_analytical_[5]));
//
////    std::cout << " \n ********************************* \n";
////    std::cout << " \n relative difference displacements [%]: \n";
////    std::vector<double> relative_err_dis(6);
////    for (int i=0; i<6; i++){
////        relative_err_dis[i]=abs(displ_analytical_[i]-displATpoints[i])/abs(displ_analytical_[i]);
////        std::cout << relative_err_dis[i] * 100 << "   ";
////    }
////    std::cout << " \n ********************************* \n";
//
//
//    // redefining the load
//    px = 155.; py = 230.; pz = 0.;
//    load_g[0] = -px; load_g[1] = -py; load_g[2] = -pz;
//    elem = 0;
//    for(int i = 0; i < mymesh.nelts * 3; i= i + 3)
//    {
//        bie::FaceData face = mesh3d.getElementData(elem);
//        il::Array2D<double> R = face.rotationMatrix(true); // g->l rotation matrix
//        load_l = il::dot(R,load_g); // get the local load
//
//        ff[i] =   load_l[0];
//        ff[i+1] = load_l[1];
//        ff[i+2] = load_l[2];
//
//        x0[i]=0.;
//        x0[i+1]=0.;
//        x0[i+2]=0.;
//
//        // incremement the element iterator
//        elem ++;
//    }
//
//    // gmres via bigwham object.
//    stat= tractionHMAT.h_IterativeSolve(ff, false, 1000, il::io, nIts, tolerance, x0);
//
//    std::cout << " Gmress success (0 ?) ";
//    if(stat == 0){std::cout << " YES "; }
//    else { std::cout << " NO ";}
//    std::cout <<"\n";
//
//    for (int i=1;i<mymesh.nelts * 3;i++ ){
//        x0_[i]=x0[i];
//    }
//
//    // get stress at obs points
//    stressATpoints = tractionHMAT.computeStresses(x0_, obsPoints, 2, properties , mymesh.coor, mymesh.conn, false);
//
//    a = 1.2; // crack radius
//    x = obsPoints[0], y = obsPoints[1], z = obsPoints[2];
//
//    stress_analytical_[0] = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    stress_analytical_[1] = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    stress_analytical_[2] = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//    stress_analytical_[3] = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py) ; // sxy
//    stress_analytical_[4] = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py) ; // sxz
//    stress_analytical_[5] = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py) ; // syz
//
//    ASSERT_NEAR(stress_analytical_[0], stressATpoints[0], abs(0.1 * stress_analytical_[0]));
//    ASSERT_NEAR(stress_analytical_[1], stressATpoints[1], abs(0.1 * stress_analytical_[1]));
//    ASSERT_NEAR(stress_analytical_[2], stressATpoints[2], abs(0.1 * stress_analytical_[2]));
//    ASSERT_NEAR(stress_analytical_[3], stressATpoints[3], abs(0.1 * stress_analytical_[3]));
//    ASSERT_NEAR(stress_analytical_[4], stressATpoints[4], abs(0.1 * stress_analytical_[4]));
//    ASSERT_NEAR(stress_analytical_[5], stressATpoints[5], abs(0.1 * stress_analytical_[5]));
//
//    x = obsPoints[3], y = obsPoints[4], z = obsPoints[5];
//    stress_analytical_[6] = bie::Sig_xx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xx_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    stress_analytical_[7] = bie::Sig_yy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_yy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    stress_analytical_[8] = bie::Sig_zz_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//    stress_analytical_[9] = bie::Sig_xy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_xy_sl_(x, y, z, a, G, nu, px, py) ; // sxy
//    stress_analytical_[10] = bie::Sig_zx_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zx_sl_(x, y, z, a, G, nu, px, py) ; // sxz
//    stress_analytical_[11] = bie::Sig_zy_nl_(x, y, z, a, G, nu, pz) + bie::Sig_zy_sl_(x, y, z, a, G, nu, px, py) ; // syz
//
//    ASSERT_NEAR(stress_analytical_[6], stressATpoints[6], abs(0.1 * stress_analytical_[6]));
//    ASSERT_NEAR(stress_analytical_[7], stressATpoints[7], abs(0.1 * stress_analytical_[7]));
//    ASSERT_NEAR(stress_analytical_[8], stressATpoints[8], abs(0.1 * stress_analytical_[8]));
//    ASSERT_NEAR(stress_analytical_[9], stressATpoints[9], abs(0.1 * stress_analytical_[9]));
//    ASSERT_NEAR(stress_analytical_[10], stressATpoints[10], abs(0.1 * stress_analytical_[10]));
//    ASSERT_NEAR(stress_analytical_[11], stressATpoints[11], abs(0.1 * stress_analytical_[11]));
//
////    std::cout << " \n ********************************* \n";
////    std::cout << " \n relative difference stresses [%]: \n";
////    for (int i=0; i<12; i++){
////        relative_err_stress[i]=abs(stress_analytical_[i]-stressATpoints[i])/abs(stress_analytical_[i]);
////        std::cout << relative_err_stress[i] * 100 << "   ";
////    }
//
//    // check displacements
//    displATpoints = tractionHMAT.computeDisplacements(x0_, obsPoints, 2, properties , mymesh.coor, mymesh.conn, false);
//
//    x = obsPoints[0], y = obsPoints[1], z = obsPoints[2];
//    displ_analytical_[0] = bie::Ux_nl_(x, y, z, a, G, nu, pz) + bie::Ux_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    displ_analytical_[1] = bie::Uy_nl_(x, y, z, a, G, nu, pz) + bie::Uy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    displ_analytical_[2] = bie::Uz_nl_(x, y, z, a, G, nu, pz) + bie::Uz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//
//    ASSERT_NEAR(displ_analytical_[0], displATpoints[0], abs(0.1 * displ_analytical_[0]));
//    ASSERT_NEAR(displ_analytical_[1], displATpoints[1], abs(0.1 * displ_analytical_[1]));
//    ASSERT_NEAR(displ_analytical_[2], displATpoints[2], abs(0.1 * displ_analytical_[2]));
//
//    x = obsPoints[3], y = obsPoints[4], z = obsPoints[5];
//    displ_analytical_[3] = bie::Ux_nl_(x, y, z, a, G, nu, pz) + bie::Ux_sl_(x, y, z, a, G, nu, px, py) ; // sxx
//    displ_analytical_[4] = bie::Uy_nl_(x, y, z, a, G, nu, pz) + bie::Uy_sl_(x, y, z, a, G, nu, px, py) ; // syy
//    displ_analytical_[5] = bie::Uz_nl_(x, y, z, a, G, nu, pz) + bie::Uz_sl_(x, y, z, a, G, nu, px, py) ; // szz
//
//    ASSERT_NEAR(displ_analytical_[3], displATpoints[3], abs(0.1 * displ_analytical_[3]));
//    ASSERT_NEAR(displ_analytical_[4], displATpoints[4], abs(0.1 * displ_analytical_[4]));
//    ASSERT_NEAR(displ_analytical_[5], displATpoints[5], abs(0.1 * displ_analytical_[5]));
//
////    std::cout << " \n ********************************* \n";
////    std::cout << " \n relative difference displacements  [%]: \n";
////    for (int i=0; i<6; i++){
////        relative_err_dis[i]=abs(displ_analytical_[i]-displATpoints[i])/abs(displ_analytical_[i]);
////        std::cout << relative_err_dis[i] * 100 << "   ";
////    }
////    std::cout << " \n ********************************* \n";
//
//}
