#ifndef BIGWHAM_BIE_ELASTIC3D_HYDROST_EIGENSTRAIN_TRIANGLE_H
#define BIGWHAM_BIE_ELASTIC3D_HYDROST_EIGENSTRAIN_TRIANGLE_H

#include <stdexcept>
#include <tuple>

#include <il/StaticArray.h>
#include <il/Array2D.h>

#include "elasticity/bie_elastostatic_eigenstrain.h"
#include "elements/boundary_element.h"

#include "elements/triangle.h"
#include "elements/rectangle.h"
#include "elements/hexahedron.hh"

#include "hydrostatic_eigenstrain_polyhedra.hh"

namespace bigwham {

// Base class full specializations required to satisfy the vtable in debug builds
// (no-op bodies; influence is always overridden by BieElastostaticEigenstrain)
template <>
std::vector<double>
BieElastostatic<Hexahedron<0>, Rectangle<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &, il::int_t, const BoundaryElement &, il::int_t) const {
    throw std::logic_error("BieElastostatic base influence called for eigenstrain kernel");
}
template <>
std::vector<double>
BieElastostatic<Tetrahedron<0>, Triangle<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &, il::int_t, const BoundaryElement &, il::int_t) const {
    throw std::logic_error("BieElastostatic base influence called for eigenstrain kernel");
}

template class BieElastostaticEigenstrain<Hexahedron<0>, Rectangle<0>, ElasticKernelType::V>;
template class BieElastostaticEigenstrain<Tetrahedron<0>, Triangle<0>, ElasticKernelType::V>;


/**
 * @brief V kernel - nuclei of strain. triangle -> segment
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
BieElastostaticEigenstrain<Hexahedron<0>, Rectangle<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    // Get the observation point
    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 3> x_obs{il::value, {
        r_col(i_r, 0), 
        r_col(i_r, 1),
        r_col(i_r, 2)
    }};

    // Get the nornal of the segment element
    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 3> n_obs{il::value, {
        n_array[0], 
        n_array[1],
        n_array[2]
    }};

    // Catsting the source_elmt to Hexahedron<0>
    const Hexahedron<0>* poly_ptr = dynamic_cast<const Hexahedron<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt could not be cast to Hexahedron<0>!");
    }

    // Compute the tractions for unit eigenstrain
    auto t_i = V_threeD_polyhedra_0(
        *poly_ptr, // Triangle<0> 
        x_obs, // const il::StaticArray<double, 3>
        n_obs, // const il::StaticArray<double, 3>
        this->elas_.shear_modulus(), 
        this->elas_.poisson_ratio()
    );

    // Return as std::vector
    std::vector<double> t_i_res = {t_i[0], t_i[1], t_i[2]};
    return t_i_res;
}

/**
 * @brief V kernel - nuclei of strain. triangle -> segment
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
BieElastostaticEigenstrain<Tetrahedron<0>, Triangle<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    // Get the observation point
    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 3> x_obs{il::value, {
        r_col(i_r, 0), 
        r_col(i_r, 1),
        r_col(i_r, 2)
    }};

    // Get the nornal of the segment element
    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 3> n_obs{il::value, {
        n_array[0], 
        n_array[1],
        n_array[2]
    }};

    // Catsting the source_elmt to Tetrahedron<0>
    const Tetrahedron<0>* poly_ptr = dynamic_cast<const Tetrahedron<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt could not be cast to Tetrahedron<0>!");
    }

    // Compute the tractions for unit eigenstrain
    auto t_i = V_threeD_polyhedra_0(
        *poly_ptr, //  Polyhedra<0> 
        x_obs, // const il::StaticArray<double, 3>
        n_obs, // const il::StaticArray<double, 3>
        this->elas_.shear_modulus(), 
        this->elas_.poisson_ratio()
    );

    // Return as std::vector
    std::vector<double> t_i_res = {t_i[0], t_i[1], t_i[2]};
    return t_i_res;
}


} // namespace bigwham

#endif // BIGWHAM_BIE_ELASTIC3D_HYDROST_EIGENSTRAIN_TRIANGLE_H
