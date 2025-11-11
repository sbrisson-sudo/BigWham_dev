#ifndef BIGWHAM_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
#define BIGWHAM_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/math.h>
#include <array>
#include <utility>

#include "elements/polygon.h"

namespace bigwham {

/**
 * @brief Integrate the plane strain eigenstrain kernel over a polygonal
 * element at a given observation point. Returns the traction given a normal.
 * 
 * @param polygon a polygon (in 2d space)
 * @param xy_obs observation point
 * @param n_obs normal
 * @param G shear modulus
 * @param nu Poisson's ratio
 * @return il::StaticArray<double, 2> tractions
 */
il::StaticArray<double, 2> V_twoD_polygon_0(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> xy_obs,
    const il::StaticArray<double, 2> n_obs,
    double G, double nu,
    bool recursive_call = false
);

/**
 * @brief Compute phi_ij = \int_V \partial_{x_i}\partial_{x_y} ln(r) dx' for 
 * a general polygon.
 * 
 * @param tri_vertices 
 * @param xy_obs 
 * @return std::array<std::array<double, 2>, 3> 
 */
std::array<std::array<double, 2>, 2>  phi_ij_polygon(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> xy_obs
);

/**
 * @brief Compute \int_A_I^A_Ip1 r_j/r^2 dx'
 * 
 * @param A_I 
 * @param A_Ip1 
 * @param xy_obs 
 * @param j 
 * @return double 
 */
std::pair<double,int> IGrad_j_lineintegral(
    const il::StaticArray<double, 2> A_I,
    const il::StaticArray<double, 2> A_Ip1,
    const il::StaticArray<double, 2> xy_obs,
    unsigned int j
);


} // namespace bigwham
#endif // BIGWHAM_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
