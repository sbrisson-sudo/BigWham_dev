
#ifndef BIGWHAM_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H
#define BIGWHAM_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/math.h>
#include <array>
#include <utility>

#include "elements/polygon.h"

namespace bigwham {

/**
 * @brief Integration method selection for axisymmetric potential computation
 */
enum class IntegrationMethod {
    QUADPACK,        // Adaptive QAGS integration (scipy quad)
    GAUSS_LEGENDRE,  // Gauss-Legendre quadrature
    FAR_FIELD        // Far-field approximation
};

/**
 * @brief Evaluate the potential phi(r) = \iiint_V dr'/|r-r'| for
 * V an axosymmetric domain of polygonal meridional cross section.
 *
 * @param polygon
 * @param rz_obs
 * @param method Integration method to use
 * @return double
 */
double phi_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    IntegrationMethod method = IntegrationMethod::QUADPACK
);

/**
 * @brief Compute and returns the 4 non-zero strain components eps_rr, eps_zz, eps_tt, eps_rz.
 *
 * @param polygon
 * @param rz_obs
 * @param nu Poisson's ratio
 * @param method Integration method to use
 * @return std::array<double, 4>
 */
std::array<double, 4> strain_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    double nu,
    IntegrationMethod method = IntegrationMethod::QUADPACK
);

/**
 * @brief Compute and returns the 4 non-zero stress components sigma_rr, sigma_zz, sigma_tt, sigma_rz.
 *
 * @param polygon
 * @param rz_obs
 * @param nu Poisson's ratio
 * @param G Shear modulus
 * @param method Integration method to use
 * @return std::array<double, 4>
 */
std::array<double, 4> stress_3daxis(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    double nu, double G,
    IntegrationMethod method = IntegrationMethod::QUADPACK
);

/**
 * @brief Returns the traction at (R,Z) projected on normal n due to the axisymmeytric
 * inclusion with polygon as its meridional cross-section.
 * 
 * @param polygon 
 * @param rz_obs 
 * @param n_obs 
 * @param G 
 * @param nu 
 * @param recursive_call 
 * @return il::StaticArray<double, 2> 
 */
il::StaticArray<double, 2> V_threeDAxis_polygon_0(
    const Polygon<0> &polygon,
    const il::StaticArray<double, 2> rz_obs,
    const il::StaticArray<double, 2> n_obs,
    double G, double nu,
    bool recursive_call = false
);

} // namespace bigwham 


#endif // BIGWHAM_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H