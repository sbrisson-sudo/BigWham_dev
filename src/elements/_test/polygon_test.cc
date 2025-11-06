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