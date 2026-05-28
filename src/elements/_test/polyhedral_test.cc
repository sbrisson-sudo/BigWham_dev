#include <gtest/gtest.h>

#include <iostream>
#include <iomanip> // For std::setw

#include <il/Array2D.h>

#include "elements/polyhedral.hh"
#include "elements/hexahedron.hh"
#include "elements/tetrahedron.hh"

using namespace bigwham;

template <typename T>
void prettyPrintArray2D(const il::Array2D<T> &M) {
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

TEST(Polyhedral, hexahedron) {

    // Unit cube
    il::Array2D<double> hex_vertices{il::value, {
        {0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0}
    }};

    Hexahedron<0> hex;

    hex.SetElement(hex_vertices);

    // std::cout << "Face normals :\n";
    // prettyPrintArray2D(hex.getFaceNormals());

    // std::cout << "Face centroids :\n";
    // prettyPrintArray2D(hex.getFaceCentroids());

    // Test geometric properties 
    ASSERT_TRUE( hex.isPointInPolyhedron({0.5, 0.5, 0.5}) );
    ASSERT_FALSE( hex.isPointInPolyhedron({1.5, 0.5, 0.5}) );

    ASSERT_FALSE( hex.isPointOnBoundary({0.5, 0.5, 0.5}) );
    ASSERT_FALSE( hex.isPointOnBoundary({1.5, 0.5, 0.5}) );
    ASSERT_TRUE( hex.isPointOnBoundary({1.0, 0.5, 0.5}) );

    double tol = hex.getTol();
    ASSERT_FALSE( hex.isPointInPolyhedron({1.0 + tol*10, 0.5, 0.5}) );
    ASSERT_TRUE( hex.isPointInPolyhedron({1.0 - tol*10, 0.5, 0.5}) );

}


TEST(Polyhedral, tetrahedron) {

    // Unit tet 
    il::Array2D<double> tet_vertices{il::value, {
        {0.0,   0.866025404,  -0.866025404,   0.0},
        {0.0,   -0.5,         -0.5,           1.0},
        {1.0,   0.0,          0.0,            0.0}
    }};

    Tetrahedron<0> tet;

    tet.SetElement(tet_vertices);

    // std::cout << "Face normals :\n";
    // prettyPrintArray2D(tet.getFaceNormals());

    // std::cout << "Face centroids :\n";
    // prettyPrintArray2D(tet.getFaceCentroids());

    // Test geometric properties 
    ASSERT_TRUE( tet.isPointInPolyhedron({0.0, 0.0, 0.25}) );
    ASSERT_FALSE( tet.isPointInPolyhedron({0.0, 0.0, 1.25}) );

    ASSERT_FALSE( tet.isPointOnBoundary({0.0, 0.0, 0.25}) );
    ASSERT_FALSE( tet.isPointOnBoundary({0.0, 0.0, 1.25}) );
    ASSERT_TRUE( tet.isPointOnBoundary({0.0, 0.0, 0.0}) );

    double tol = tet.getTol();
    ASSERT_TRUE( tet.isPointInPolyhedron({0.0, 0.0, 0.0 + tol*10}) );
    ASSERT_FALSE( tet.isPointInPolyhedron({0.0, 0.0, 0.0 - tol*10}) );

}