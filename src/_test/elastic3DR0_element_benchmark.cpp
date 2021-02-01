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
#include "elastic3DR0_element_benchmark.h"
#include <BigWham.h>

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
    const std::string tractionKernel = "3DR0_traction";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points
    std::vector<double> dd_solution(64 * 3), obsPoints = {1., 2., 3., -1., -2., -3.};

    for (int x = 0; x < 64 * 3; ++x) dd_solution[x] = 1.;

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

    for (int x = 0; x < 64 * 3; ++x) dd_solution[x] = 1.;

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
    const std::string tractionKernel = "3DR0_traction";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points
    std::vector<double> dd_solution(64 * 3), obsPoints = {1., 2., 3., -1.,-2.,-3.};

    for (int x = 0; x < 64 * 3; ++x) dd_solution[x] = 1.;

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

    for (int x = 0; x < 64 * 3; ++x) dd_solution[x] = 1.;

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
     * Then we apply a rotation to te plane according to the following matrix:
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
    const std::string tractionKernel = "3DR0_traction";

    Bigwhamio tractionHMAT;
    tractionHMAT.set(mymesh.coor, mymesh.conn, tractionKernel, properties, max_leaf_size, eta, eps_aca);

    // get stress at points by prescribing global DD
    std::vector<double>  dd_solution_g(64 * 3), obsPoints = {0.,0.,0., 100., 100., 100.};

    for(int x = 0; x < 64*3; x=x+3)
    { dd_solution_g[x] = -0.5687111245916712;
      dd_solution_g[x+1] = -0.09950371902099886;
      dd_solution_g[x+2] = -1.6329931618554523;}

    std::vector<double> stressATpoints_0 = tractionHMAT.computeStresses(dd_solution_g, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, true);

    // get stress at points by prescribing local DD
    std::vector<double>  dd_solution_l(64 * 3);

    for(int x = 0; x < 64*3; x=x+3)
    { dd_solution_l[x] = 1;
      dd_solution_l[x+1] = 1;
      dd_solution_l[x+2] = 1;}

    std::vector<double> stressATpoints = tractionHMAT.computeStresses(dd_solution_l, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, false);

    // check that, if properly done, prescribing global DD or local DDs does not matter
    ASSERT_NEAR(stressATpoints_0[0], stressATpoints[0], 2.e-5);
    ASSERT_NEAR(stressATpoints_0[1], stressATpoints[1], 2.e-5);
    ASSERT_NEAR(stressATpoints_0[2], stressATpoints[2], 2.e-5);
    ASSERT_NEAR(stressATpoints_0[3], stressATpoints[3], 2.e-5);
    ASSERT_NEAR(stressATpoints_0[4], stressATpoints[4], 2.e-5);
    ASSERT_NEAR(stressATpoints_0[5], stressATpoints[5], 2.e-5);


    // solution at point {0,0,0}:

    // Rotation matricies:
    il::Array2D<double> R{3,3,0.}; // global to local ...R(g->l)
    R(0,0) = -0.690578; R(0,1) = 0.597022; R(0,2) = -0.408248;
    R(1,0) = 0.568711;  R(1,1) = 0.0995037; R(1,2) = -0.816497;
    R(2,0) = -0.446844; R(2,1) = -0.79603; R(2,2) = -0.408248;

    il::Array2D<double> RT{3,3,0.}; // local to global ...R(l->g)
    RT(0,0) = -0.690578; RT(0,1) = 0.568711;  RT(0,2) = -0.446844;
    RT(1,0) =  0.597022; RT(1,1) = 0.0995037; RT(1,2) = -0.79603;
    RT(2,0) =  -0.408248; RT(2,1) = -0.816497;  RT(2,2) =  -0.408248;

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
    for(int x = 0; x < 64*3; x=x+3)
    {   dd_solution_g[x] = -0.5687111245916712;
        dd_solution_g[x+1] = -0.09950371902099886;
        dd_solution_g[x+2] = -1.6329931618554523;}

    std::vector<double> displATpoints_0 = displacementHMAT.computeDisplacements(dd_solution_g, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, true);

    // get stress at points by prescribing local DD
    for(int x = 0; x < 64*3; x=x+3)
    { dd_solution_l[x] = 1;
        dd_solution_l[x+1] = 1;
        dd_solution_l[x+2] = 1;}

    std::vector<double> displATpoints = displacementHMAT.computeDisplacements(dd_solution_l, obsPoints, 2, properties ,mymesh.coor, mymesh.conn, false);

    // check that, if properly done, prescribing global DD or local DDs does not matter
    ASSERT_NEAR(displATpoints_0[0], displATpoints[0], 2.e-5);
    ASSERT_NEAR(displATpoints_0[1], displATpoints[1], 2.e-5);
    ASSERT_NEAR(displATpoints_0[2], displATpoints[2], 2.e-5);
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


