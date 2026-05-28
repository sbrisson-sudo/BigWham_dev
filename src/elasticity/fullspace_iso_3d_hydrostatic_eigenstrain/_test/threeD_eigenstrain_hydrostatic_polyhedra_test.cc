#include <gtest/gtest.h>

#include <iostream>
#include <iomanip> // For std::setw

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/Array2D.h>

// #include "elements/triangle.h"
#include "elements/hexahedron.hh"
#include "elements/tetrahedron.hh"

#include "elasticity/fullspace_iso_3d_hydrostatic_eigenstrain/hydrostatic_eigenstrain_polyhedra.hh"

using namespace bigwham;

template <typename T>
void prettyPrintArray(const il::Array2D<T> &M) {
    // Get the dimensions of the matrix
    il::int_t rows = M.size(0);
    il::int_t cols = M.size(1);

    // Determine the maximum width needed for any element in the matrix
    int maxWidth = 0;
    for (il::int_t i = 0; i < rows; ++i) {
        for (il::int_t j = 0; j < cols; ++j) {
            std::ostringstream oss;
            oss << M(i, j);
            int width = oss.str().length();
            if (width > maxWidth) {
                maxWidth = width;
            }
        }
    }

    // Print the matrix with proper formatting
    for (il::int_t i = 0; i < rows; ++i) {
        for (il::int_t j = 0; j < cols; ++j) {
            std::cout << std::setw(maxWidth + 2) << M(i, j);
        }
        std::cout << std::endl;
    }
}

template <typename T, long int p>
void prettyPrintArray(const il::StaticArray<T, p> &v) {

    // Determine the maximum width needed for any element in the matrix
    int maxWidth = 0;
    for (il::int_t i = 0; i < p; ++i) {
        std::ostringstream oss;
        oss << v[i];
        int width = oss.str().length();
        if (width > maxWidth) {
            maxWidth = width;
        }
    }

    // Print the vector with proper formatting
    for (il::int_t i = 0; i < p; ++i) 
        std::cout << std::setw(maxWidth + 2) << v[i];
        
    std::cout << std::endl;
    
}

TEST(ThreeDHex0R0_V, test1) {

    double tol = 0.01;

    double nu = 0.25;
    double G = 1.0;

    // Unit cube
    il::Array2D<double> hex_vertices{il::value, {
        {0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0}
    }};

    Hexahedron<0> hex;

    hex.SetElement(hex_vertices);

    il::StaticArray<double, 3> x_obs{0.};
    il::StaticArray<double, 3> n_obs{il::value, {1.0, 0., 0.}};

    // In the cube 
    x_obs[0] = 0.5;
    x_obs[1] = 0.5;
    x_obs[2] = 0.5;

    // std::cout << "--------------------\n";
    // std::cout << "x_obs = ";
    // prettyPrintArray(x_obs);

    auto t_i = V_threeD_polyhedra_0(hex, x_obs, n_obs, G, nu);

    std::vector<double> t_i_ref_1{-2.22222222, 0., 0.};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_1[i], tol);
    }

    // outside the cube 
    x_obs[0] = 1.5;
    x_obs[1] = 0.5;
    x_obs[2] = 0.5;

    // std::cout << "--------------------\n";
    // std::cout << "x_obs = ";
    // prettyPrintArray(x_obs);

    t_i = V_threeD_polyhedra_0(hex, x_obs, n_obs, G, nu);

    std::vector<double> t_i_ref_2{-0.44927462, 0., 0.};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_2[i], tol);
    }

    // on the boundary of the cube 
    x_obs[0] = 1.0;
    x_obs[1] = 0.5;
    x_obs[2] = 0.5;

    // std::cout << "--------------------\n";
    // std::cout << "x_obs = ";
    // prettyPrintArray(x_obs);

    t_i = V_threeD_polyhedra_0(hex, x_obs, n_obs, G, nu);

    std::vector<double> t_i_ref_3{-1.45301927, 0., 0.};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_3[i], tol);
    }
}

TEST(ThreeDTet0T0_V, test1) {

    double tol = 1e-3;

    double nu = 0.25;
    double G = 1.0;

    // Regular tetrahedron-like shape
    // v0=(0,0,h), v1=(cos(pi/6), -sin(pi/6), 0), v2=(-cos(pi/6), -sin(pi/6), 0), v3=(0,1,0)
    // h=1, cos(pi/6)=sqrt(3)/2, sin(pi/6)=0.5
    const double c = std::sqrt(3.0) / 2.0; // cos(pi/6)
    const double s = 0.5;                   // sin(pi/6)

    il::Array2D<double> tet_vertices{il::value, {
        {0.0,  c, -c, 0.0},  // x coordinates of v0,v1,v2,v3
        {0.0, -s, -s, 1.0},  // y coordinates
        {1.0,  0.0,  0.0, 0.0}   // z coordinates
    }};

    Tetrahedron<0> tet;
    tet.SetElement(tet_vertices);

    il::StaticArray<double, 3> n_obs{il::value, {0.0, 0.0, 1.0}};
    il::StaticArray<double, 3> x_obs{0.};

    // Inside
    x_obs[0] = 0.0; x_obs[1] = 0.0; x_obs[2] = 0.5;
    auto t_i = V_threeD_polyhedra_0(tet, x_obs, n_obs, G, nu);
    std::vector<double> t_i_ref_1{0.0, 0.0, -1.77676398};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_1[i], tol);
    }

    // On boundary
    x_obs[0] = 0.0; x_obs[1] = 0.0; x_obs[2] = 0.0;
    t_i = V_threeD_polyhedra_0(tet, x_obs, n_obs, G, nu);
    std::vector<double> t_i_ref_2{0.0, 0.0, -1.48453059};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_2[i], tol);
    }

    // Outside
    x_obs[0] = 0.0; x_obs[1] = 0.0; x_obs[2] = -0.5;
    t_i = V_threeD_polyhedra_0(tet, x_obs, n_obs, G, nu);
    std::vector<double> t_i_ref_3{0.0, 0.0, -0.67974900};
    for (int i(0); i<3; i++){
        ASSERT_NEAR(t_i[i], t_i_ref_3[i], tol);
    }
}
