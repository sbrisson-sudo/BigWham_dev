//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.05.2026.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#ifndef BIGWHAM_BIE_ELASTOSTATIC_AXI3D_MODE1_H
#define BIGWHAM_BIE_ELASTOSTATIC_AXI3D_MODE1_H

#include "elasticity/bie_elastostatic.h"

namespace bigwham {

// Mode-I only axisymmetric ring kernel.
// DOF dimension = 1 (normal/opening displacement discontinuity only).
// Analogous to BieElastostaticModeI but for the axisymmetric ring geometry.
template <class Es, class Er, ElasticKernelType k>
class BieElastostaticAxi3DModeI : public BieElastostatic<Es, Er, k> {
    using BieElastostatic<Es, Er, k>::BieElastostatic;

public:
    BieElastostaticAxi3DModeI() : BieElastostatic<Es, Er, k>() {};

    BieElastostaticAxi3DModeI(bigwham::ElasticProperties &elas, il::int_t dim)
        : BieElastostatic<Es, Er, k>() {
        IL_EXPECT_FAST(dim == 2);
        this->elas_ = elas;
        this->dof_dimension_ = 1;
        this->spatial_dimension_ = dim;
    };

    BieElastostaticAxi3DModeI(bigwham::ElasticProperties &elas, il::int_t dim,
                               bool local_unknowns, bool local_co_variables)
        : BieElastostatic<Es, Er, k>() {
        IL_EXPECT_FAST(dim == 2);
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
#endif // BIGWHAM_BIE_ELASTOSTATIC_AXI3D_MODE1_H
