#ifndef BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
#define BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H

#include <tuple>

#include <il/StaticArray.h>
#include <il/Array2D.h>

#include "elasticity/bie_elastostatic_eigenstrain.h"
#include "elements/boundary_element.h"
#include "elements/segment.h"
#include "elements/triangle.h"
#include "elastic_2dT0S0_V_element.hh"

namespace bigwham {

template class BieElastostaticEigenstrain<Triangle<0>, Segment<0>, ElasticKernelType::V>;


//  V kernel - nuclei of strain
template <>
std::vector<double>
BieElastostaticEigenstrain<Triangle<0>, Segment<0>, ElasticKernelType::V>::influence(
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

    // Get the triangle vertices 
    auto tri_vertices_array = source_elt.vertices();
    std::array<std::array<double, 2>, 3> tri_vertices{{
        {tri_vertices_array(0,0), tri_vertices_array(0,1)},
        {tri_vertices_array(1,0), tri_vertices_array(1,1)},
        {tri_vertices_array(2,0), tri_vertices_array(2,1)}
    }};

    // Compute the tractions for unit eigenstrain
    auto t_i = V_twoD_triangle_0(
        tri_vertices, // const std::array<std::array<double, 2>, 3>  
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

#endif // BIGWHAM_BIE_ELASTIC2D_HYDROST_EIGENSTRAIN_TRIANGLE_H
