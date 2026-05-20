//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 14.05.2026.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#ifndef BIGWHAM_BIE_ELASTOSTATIC_MODE1_SYM_H
#define BIGWHAM_BIE_ELASTOSTATIC_MODE1_SYM_H

#include "elasticity/bie_elastostatic.h"

namespace bigwham {

// Symmetric version of BieElastostaticModeI with symmetry plane x=0 (y-axis).
// The mesh is defined on x>0; the influence function sums contributions from
// both the actual element and its mirror image reflected about x=0.
template <class Es, class Er, ElasticKernelType k>
class BieElastostaticModeISym : public BieElastostatic<Es, Er, k> {
    using BieElastostatic<Es, Er, k>::BieElastostatic;

public:
    BieElastostaticModeISym() : BieElastostatic<Es, Er, k>() {};
    BieElastostaticModeISym(bigwham::ElasticProperties &elas, il::int_t dim) :
    BieElastostatic<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = 1;
        this->spatial_dimension_ = dim;
    };

    BieElastostaticModeISym(bigwham::ElasticProperties &elas,
        il::int_t dim, bool local_unknowns, bool local_co_variables) :
    BieElastostatic<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = 1;
        this->spatial_dimension_ = dim;
        this->local_unknowns_ = local_unknowns;
        this->local_co_variables_ = local_co_variables;
    };

    std::vector<double> influence(const BoundaryElement &, il::int_t,
                                  const BoundaryElement &,
                                  il::int_t) const override;
};

} // namespace bigwham

#endif // BIGWHAM_BIE_ELASTOSTATIC_MODE1_SYM_H
