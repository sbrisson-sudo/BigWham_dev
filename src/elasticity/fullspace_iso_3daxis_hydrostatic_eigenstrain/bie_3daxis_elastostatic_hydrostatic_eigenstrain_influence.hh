#ifndef BIGWHAM_BIE_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H
#define BIGWHAM_BIE_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H

#include <tuple>

#include <il/StaticArray.h>
#include <il/Array2D.h>

#include "elasticity/bie_elastostatic_eigenstrain.h"
#include "elements/boundary_element.h"
#include "elements/segment.h"
#include "elements/triangle.h"
#include "elements/rectangle.h"
#include "hydrostatic_eigenstrain_3daxis.hh"

namespace bigwham {

template class BieElastostaticEigenstrainAxi3D<Triangle<0>, Segment<0>, ElasticKernelType::V>;
template class BieElastostaticEigenstrainAxi3D<Rectangle<0>, Segment<0>, ElasticKernelType::V>;


/**
 * @brief V kernel - nuclei of strain. triangle axis. -> segment
 * 
 * @tparam  
 * @param source_elt 
 * @param i_s 
 * @param receiver_elt 
 * @param i_r 
 * @return std::vector<double> 
 */
template <>
inline std::vector<double>
BieElastostaticEigenstrainAxi3D<Triangle<0>, Segment<0>, ElasticKernelType::V>::influence(
        const BoundaryElement &source_elt, il::int_t i_s,
        const BoundaryElement &receiver_elt, il::int_t i_r) const {

    // Get the observation point
    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 2> xy_obs{il::value, {
        r_col(i_r, 0), 
        r_col(i_r, 1)}
    };

    // Get the nornal of the segment element
    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 2> n_obs{il::value, {
        n_array[0], 
        n_array[1]}
    };

    // Catsting the source_elmt to Triangle<0>
    const Triangle<0>* poly_ptr = dynamic_cast<const Triangle<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt is not a Triangle<0>!");
    }

    // Compute the tractions for unit eigenstrain
    auto t_i = V_threeDAxis_polygon_0(
        *poly_ptr, // Triangle<0> 
        xy_obs, // const il::StaticArray<double, 2>
        n_obs, // const il::StaticArray<double, 2>
        this->elas_.shear_modulus(), 
        this->elas_.poisson_ratio()
    );

    // Return as std::vector
    std::vector<double> t_i_res = {t_i[0], t_i[1]};
    return t_i_res;
}

/**
 * @brief V kernel - nuclei of strain. rectangle -> segment
 * 
 * @tparam  
 * @param source_elt 
 * @param i_s 
 * @param receiver_elt 
 * @param i_r 
 * @return std::vector<double> 
 */
template <>
std::vector<double>
BieElastostaticEigenstrainAxi3D<Rectangle<0>, Segment<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    // Get the observation point
    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 2> xy_obs{il::value, {
        r_col(i_r, 0), 
        r_col(i_r, 1)}
    };

    // Get the nornal of the segment element
    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 2> n_obs{il::value, {
        n_array[0], 
        n_array[1]}
    };

    // Catsting the source_elmt to Triangle<0>
    const Rectangle<0>* poly_ptr = dynamic_cast<const Rectangle<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt is not a Rectangle<0>!");
    }

    // Compute the tractions for unit eigenstrain
    auto t_i = V_threeDAxis_polygon_0(
        *poly_ptr, // Rectangle<0> 
        xy_obs, // const il::StaticArray<double, 2>
        n_obs, // const il::StaticArray<double, 2>
        this->elas_.shear_modulus(), 
        this->elas_.poisson_ratio()
    );

    // Return as std::vector
    std::vector<double> t_i_res = {t_i[0], t_i[1]};
    return t_i_res;
}


} // namespace bigwham

#endif // BIGWHAM_BIE_ELASTIC3DAXIS_HYDROST_EIGENSTRAIN_H
