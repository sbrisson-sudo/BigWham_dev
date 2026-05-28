#ifndef BIGWHAM_ELASTIC_HYDROSTATIC_EIGENSTRAIN_POLYHEDRA
#define BIGWHAM_ELASTIC_HYDROSTATIC_EIGENSTRAIN_POLYHEDRA

#include <array>

#include <il/StaticArray.h>
#include <il/Array2D.h>
#include <il/StaticArray2D.h>
#include "elements/polyhedral.hh"


namespace bigwham {

/**
 * @brief Integrate the plane strain eigenstrain kernel over a polygonal
 * element at a given observation point. Returns the traction given a normal.
 * 
 * @param polyhedra a polyhedra
 * @param xy_obs observation point
 * @param n_obs normal
 * @param G shear modulus
 * @param nu Poisson's ratio
 * @return il::StaticArray<double, 2> tractions
 */
il::StaticArray<double, 3> V_threeD_polyhedra_0(
    const Polyhedral<0> &polyhedra,
    const il::StaticArray<double, 3> xy_obs,
    const il::StaticArray<double, 3> n_obs,
    double G, double nu,
    bool recursive_call = false
);

/**
 * @brief Compute phi_ij = \int_V \partial_{x_i}\partial_{x_y} 1/r dx' for 
 * a general polyhedra.
 * 
 * @param polyhedra 
 * @param xy_obs 
 * @return std::array<std::array<double, 2>, 3> 
 */
std::array<std::array<double, 3>, 3>  phi_ij_polyhedra(
    const Polyhedral<0> &polyhedra,
    const il::StaticArray<double, 3> xy_obs
);



} // namespace bigwham






#endif // BIGWHAM_ELASTIC_HYDROSTATIC_EIGENSTRAIN_POLYHEDRA