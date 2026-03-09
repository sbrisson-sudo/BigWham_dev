//
// This file is part of BigWham.
//
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_BIEELASTOSTATIC_EIGENSTRAIN_H
#define BIGWHAM_BIEELASTOSTATIC_EIGENSTRAIN_H

#include <vector>

#include "elasticity/bie_elastostatic.h"

namespace bigwham {

template <class Es, class Er,ElasticKernelType k>
class BieElastostaticEigenstrain : public BieElastostatic<Es, Er, k> {
    using BieElastostatic<Es, Er, k>::BieElastostatic;

public:
    BieElastostaticEigenstrain() : BieElastostatic<Es, Er, k>() {};
    BieElastostaticEigenstrain(bigwham::ElasticProperties& elas,il::int_t dim) :
    BieElastostatic<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = {il::value, {dim, 1}};
        this->spatial_dimension_ = dim;
    };

    BieElastostaticEigenstrain(bigwham::ElasticProperties&elas,
    il::int_t dim, bool local_unknowns, bool local_co_variables):BieElastostatic<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = {il::value, {dim, 1}};
        this->spatial_dimension_ = dim;
        this->local_unknowns_ = local_unknowns;
        this->local_co_variables_ = local_co_variables;
    };

    std::vector<double> influence(const BoundaryElement &, il::int_t,
                                  const BoundaryElement &,
                                  il::int_t) const override;

};

// a dummy derived class for axi P0 kernel....
template <class Es, class Er, ElasticKernelType k>
class BieElastostaticEigenstrainAxi3D : public BieElastostaticEigenstrain<Es, Er, k> {
    using BieElastostaticEigenstrain<Es, Er, k>::BieElastostaticEigenstrain;

public:
    BieElastostaticEigenstrainAxi3D() : BieElastostaticEigenstrain<Es, Er, k>() {};
    BieElastostaticEigenstrainAxi3D(bigwham::ElasticProperties& elas,il::int_t dim) :
    BieElastostaticEigenstrain<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = {il::value, {dim, 1}};
        this->spatial_dimension_ = dim;
    };

    BieElastostaticEigenstrainAxi3D(bigwham::ElasticProperties&elas,
    il::int_t dim, bool local_unknowns, bool local_co_variables):BieElastostaticEigenstrain<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = {il::value, {dim, 1}};
        this->spatial_dimension_ = dim;
        this->local_unknowns_ = local_unknowns;
        this->local_co_variables_ = local_co_variables;
    };

    std::vector<double> influence(const BoundaryElement &, il::int_t,
                                  const BoundaryElement &,
                                  il::int_t) const override;
};

} // namespace bigwham

#endif // BIGWHAM_BIEELASTOSTATIC_EIGENSTRAIN_H
