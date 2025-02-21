//
// This file is part of BigWham.
//
// Created by Alexis Sáez Uribe on 30/05/2022.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne),
// Switzerland, Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//
//
// Edited by Ankit on 16 March 2023
// modified new API by Brice on 4 July 2024

#ifndef BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
#define BIGWHAM_ELASTICAXI3DP0_ELEMENT_H

#include <vector>

#include <il/StaticArray2D.h>

#include "core/elastic_properties.h"
#include "elements/boundary_element.h"
#include "elasticity/bie_elastostatic.h"
#include "elements/segment.h"
#include "elastic_axi3dP0_element.h"

namespace bigwham {

// a dummy derived class for axi P0 kernel....
template <class Es, class Er, ElasticKernelType k>
class BieElastostaticAxi3D : public BieElastostatic<Es, Er, k> {
    using BieElastostatic<Es, Er, k>::BieElastostatic;

public:
    BieElastostaticAxi3D() : BieElastostatic<Es, Er, k>(){};

    BieElastostaticAxi3D(bigwham::ElasticProperties &elas, il::int_t dim)
        : BieElastostatic<Es, Er, k>() {
            IL_EXPECT_FAST(dim == 2);
            this->elas_ = elas;
            this->dof_dimension_ = dim;
            this->spatial_dimension_ = dim;
            this->kernel_properties_ = il::Array<double>(1, 10000.0);
    };

    BieElastostaticAxi3D(bigwham::ElasticProperties &elas, il::int_t dim,bool local_unknowns, bool local_co_variables)
                : BieElastostatic<Es, Er, k>() {
        IL_EXPECT_FAST(dim == 2);
        this->elas_ = elas;
        this->dof_dimension_ = dim;
        this->spatial_dimension_ = dim;
        this->kernel_properties_ = il::Array<double>(1, 10000.0);
        this->local_unknowns_ = local_unknowns;
        this->local_co_variables_ = local_co_variables;
    };

    std::vector<double> influence(const BoundaryElement &, il::int_t,
                                  const BoundaryElement &,
                                  il::int_t) const override;
};

//////////////////////////////////////////////////////////////////////////////////////////////
template <>
inline std::vector<double>
BieElastostaticAxi3D<Segment<0>, Segment<0>, ElasticKernelType::H>::influence(
        const BoundaryElement &source_elt, il::int_t i_s,
        const BoundaryElement &receiver_elt, il::int_t i_r) const {
    // get constitutive parameters. NOTE:   Poisson's ratio is taken as zero
    double G = this->elas_.shear_modulus();

    auto src_vertices = source_elt.vertices();

    double rExt = std::sqrt(src_vertices(1, 0) * src_vertices(1, 0) +
                            src_vertices(1, 1) * src_vertices(1, 1));
    double rInt = std::sqrt(src_vertices(0, 0) * src_vertices(0, 0) +
                            src_vertices(0, 1) * src_vertices(0, 1));

    // std::cout << "Vertices good \n";

    // double rExt = source_elt.Xmid(0) + source_elt.size() / 2.0;
    // double rInt = source_elt.Xmid(0) - source_elt.size() / 2.0;

    auto rec_colpts = receiver_elt.collocation_points();
    double rObs = std::sqrt(rec_colpts(0, 0) * rec_colpts(0, 0) +
                            rec_colpts(0, 1) * rec_colpts(0, 1));

    // std::cout << "Collocation good \n";

    // double rObs = receiver_elt.Xmid(0);

    double IF = stress_disk_dislocation(rObs, rExt) - stress_disk_dislocation(rObs, rInt);

    il::StaticArray2D<double, 2, 2> traction_vector; // (shear, normal)

    // shear stress
    //  effect of shear dd
    traction_vector(0, 0) = 2.0 * G * IF;
    //  effect of normal dd
    traction_vector(0, 1) = 0.0;

    // normal stress
    //  effect of shear dd
    traction_vector(1, 0) = 0.0;
    //  effect of normal dd
    traction_vector(1, 1) = 2.0 * G * IF; // 2G=E when nu=0

    // std vector output in column major format
    std::vector<double> stnl(4, 0.);
    int k = 0;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            stnl[k] = traction_vector(i, j);
            k++;
        }
    }
    return stnl;
}



} // namespace bigwham

#endif // BIGWHAM_ELASTICAXI3DP0_ELEMENT_H
