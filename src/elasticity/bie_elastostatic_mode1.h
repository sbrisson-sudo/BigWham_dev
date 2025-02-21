//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.01.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//

#pragma once
#ifndef BIGWHAM_BIE_ELASTOSTATIC_MODE1_H
#define BIGWHAM_BIE_ELASTOSTATIC_MODE1_H

#include "elasticity/bie_elastostatic.h"

namespace bigwham {

template <class Es, class Er,ElasticKernelType k>
class BieElastostaticModeI : public BieElastostatic<Es, Er, k> {
    using BieElastostatic<Es, Er, k>::BieElastostatic;

public:
    BieElastostaticModeI() : BieElastostatic<Es, Er, k>() {};
    BieElastostaticModeI(bigwham::ElasticProperties& elas,il::int_t dim) :
    BieElastostatic<Es, Er, k>() {
        this->elas_ = elas;
        this->dof_dimension_ = 1;
        this->spatial_dimension_ = dim;
    };

    BieElastostaticModeI(bigwham::ElasticProperties&elas,
    il::int_t dim, bool local_unknowns, bool local_co_variables):BieElastostatic<Es, Er, k>() {
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
}
#endif //BIGWHAM_BIE_ELASTOSTATIC_MODE1_H
