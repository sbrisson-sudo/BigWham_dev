#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include "elements/polygon.h"
#include "elements/triangle.h"
#include "elements/rectangle.h"

TEST(Polygon, triangle_point_geom_relations) {

    // Column major !!
    il::Array2D<double> tri_vertices{il::value, {
        {1.0,       -0.5,       -0.5},
        {0.0,       0.8660254,  -0.8660254},
        {0.0,       0.0,        0.0}
    }};    
    
    bigwham::Triangle<0> triangle;
    triangle.SetElement(tri_vertices);

    ASSERT_TRUE(triangle.isPointInPolygon({0., 0.})); // In polygon
    ASSERT_FALSE(triangle.isPointInPolygon({1.5, 0.})); // Outside polygon

    ASSERT_TRUE(triangle.isPointOnBoundary({-0.5, 0.})); // On Boundary
    ASSERT_FALSE(triangle.isPointOnBoundary({0., 0.})); // Not on boundary
    ASSERT_FALSE(triangle.isPointOnBoundary({1.5, 0.})); // Not on boundary
}

TEST(Polygon, rectangle_point_geom_relations) {

    // Column major !!
    il::Array2D<double> rec_vertices{il::value, {
        {0.0,   1.0,    1.0,    0.0},
        {0.0,   0.0,    1.0,    1.0},
        {0.0,   0.0,    0.0,    0.0}
    }};    
    
    bigwham::Rectangle<0> rectangle;
    rectangle.SetElement(rec_vertices);

    ASSERT_TRUE(rectangle.isPointInPolygon({0.5, 0.5})); // In polygon
    ASSERT_FALSE(rectangle.isPointInPolygon({1.5, 0.})); // Outside polygon

    ASSERT_TRUE(rectangle.isPointOnBoundary({0.5, 0.})); // On Boundary
    ASSERT_FALSE(rectangle.isPointOnBoundary({0.5, 0.5})); // Not on boundary
    ASSERT_FALSE(rectangle.isPointOnBoundary({1.5, 0.})); // Not on boundary

    // Test tolerance 
    double tol = rectangle.getTol();
    ASSERT_TRUE(rectangle.isPointInPolygon({1.0-2*tol, 0.5})); // In polygon
    ASSERT_TRUE(rectangle.isPointInPolygon({1.0+0.5*tol, 0.5})); // Outside polygon but within tol
    ASSERT_FALSE(rectangle.isPointInPolygon({1.0+2*tol, 0.5})); // Outside polygon
}

TEST(Polygon, rectangle_area) {
    // Test area calculation for a rectangle using shoelace formula
    // Rectangle with width = 1.0 and height = 1.0
    il::Array2D<double> rec_vertices{il::value, {
        {0.0,   1.0,    1.0,    0.0},
        {0.0,   0.0,    1.0,    1.0},
        {0.0,   0.0,    0.0,    0.0}
    }};

    bigwham::Rectangle<0> rectangle;
    rectangle.SetElement(rec_vertices);

    // Area should be 1.0 * 1.0 = 1.0
    ASSERT_NEAR(rectangle.size(), 1.0, 1e-12);

    // Test with a different size rectangle: 2.5 x 3.7
    il::Array2D<double> rec_vertices2{il::value, {
        {0.0,   2.5,    2.5,    0.0},
        {0.0,   0.0,    3.7,    3.7},
        {0.0,   0.0,    0.0,    0.0}
    }};

    bigwham::Rectangle<0> rectangle2;
    rectangle2.SetElement(rec_vertices2);

    // Area should be 2.5 * 3.7 = 9.25
    ASSERT_NEAR(rectangle2.size(), 9.25, 1e-12);
}

TEST(Polygon, triangle_area) {
    // Test area calculation for a triangle using shoelace formula
    // Right triangle with base = 1.0 and height = 0.8660254 (sqrt(3)/2)
    il::Array2D<double> tri_vertices{il::value, {
        {1.0,       -0.5,       -0.5},
        {0.0,       0.8660254,  -0.8660254},
        {0.0,       0.0,        0.0}
    }};

    bigwham::Triangle<0> triangle;
    triangle.SetElement(tri_vertices);

    // For this equilateral triangle with side length 1.0:
    // Area = sqrt(3)/4 * side^2 = sqrt(3)/4 * (1.5)^2 ≈ 0.974278...
    // Actually, let's compute it exactly from the vertices
    // Using shoelace: A = 0.5 * |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
    // = 0.5 * |1.0*(0.8660254-(-0.8660254)) + (-0.5)*(-0.8660254-0.0) + (-0.5)*(0.0-0.8660254)|
    // = 0.5 * |1.0*1.7320508 + (-0.5)*(-0.8660254) + (-0.5)*(-0.8660254)|
    // = 0.5 * |1.7320508 + 0.4330127 + 0.4330127|
    // = 0.5 * 2.5980762 = 1.2990381
    double expected_area = 1.2990381;
    ASSERT_NEAR(triangle.size(), expected_area, 1e-6);

    // Test with a simple right triangle at origin
    // Vertices: (0,0), (1,0), (0,1)
    il::Array2D<double> tri_vertices2{il::value, {
        {0.0,   1.0,    0.0},
        {0.0,   0.0,    1.0},
        {0.0,   0.0,    0.0}
    }};

    bigwham::Triangle<0> triangle2;
    triangle2.SetElement(tri_vertices2);

    // Area should be 0.5 * base * height = 0.5 * 1.0 * 1.0 = 0.5
    ASSERT_NEAR(triangle2.size(), 0.5, 1e-12);
}