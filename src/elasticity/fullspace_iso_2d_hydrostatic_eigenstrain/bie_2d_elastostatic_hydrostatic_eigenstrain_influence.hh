#ifndef BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
#define BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H

#include <stdexcept>
#include <tuple>

#include <il/StaticArray.h>
#include <il/Array2D.h>

#include "elasticity/bie_elastostatic_eigenstrain.h"
#include "elements/boundary_element.h"
#include "elements/segment.h"
#include "elements/triangle.h"
#include "elements/rectangle.h"
#include "hydrostatic_eigenstrain_polygon.hh"

namespace bigwham {

// Base class full specializations required to satisfy the vtable in debug builds
// (no-op bodies; influence is always overridden by BieElastostaticEigenstrain)
template <>
std::vector<double>
BieElastostatic<Triangle2D<0>, Segment<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &, il::int_t, const BoundaryElement &, il::int_t) const {
    throw std::logic_error("BieElastostatic base influence called for eigenstrain kernel");
}
template <>
std::vector<double>
BieElastostatic<Rectangle2D<0>, Segment<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &, il::int_t, const BoundaryElement &, il::int_t) const {
    throw std::logic_error("BieElastostatic base influence called for eigenstrain kernel");
}

template class BieElastostaticEigenstrain<Triangle2D<0>, Segment<0>, ElasticKernelType::V>;
template class BieElastostaticEigenstrain<Rectangle2D<0>, Segment<0>, ElasticKernelType::V>;


/**
 * @brief V kernel - nuclei of strain. triangle (2D) -> segment
 */
template <>
std::vector<double>
BieElastostaticEigenstrain<Triangle2D<0>, Segment<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 2> xy_obs{il::value, {r_col(i_r, 0), r_col(i_r, 1)}};

    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 2> n_obs{il::value, {n_array[0], n_array[1]}};

    const Triangle2D<0>* poly_ptr = dynamic_cast<const Triangle2D<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt is not a Triangle2D<0>!");
    }

    auto t_i = V_twoD_polygon_0(
        *poly_ptr,
        xy_obs,
        n_obs,
        this->elas_.shear_modulus(),
        this->elas_.poisson_ratio()
    );

    return {t_i[0], t_i[1]};
}

/**
 * @brief V kernel - nuclei of strain. rectangle (2D) -> segment
 */
template <>
std::vector<double>
BieElastostaticEigenstrain<Rectangle2D<0>, Segment<0>, ElasticKernelType::V>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    auto r_col = receiver_elt.collocation_points();
    il::StaticArray<double, 2> xy_obs{il::value, {r_col(i_r, 0), r_col(i_r, 1)}};

    auto n_array = receiver_elt.normal();
    il::StaticArray<double, 2> n_obs{il::value, {n_array[0], n_array[1]}};

    const Rectangle2D<0>* poly_ptr = dynamic_cast<const Rectangle2D<0>*>(&source_elt);
    if (!poly_ptr) {
        throw std::runtime_error("Error: source_elt is not a Rectangle2D<0>!");
    }

    auto t_i = V_twoD_polygon_0(
        *poly_ptr,
        xy_obs,
        n_obs,
        this->elas_.shear_modulus(),
        this->elas_.poisson_ratio()
    );

    return {t_i[0], t_i[1]};
}


} // namespace bigwham

#endif // BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
