//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.05.2026.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once

#include <vector>

#include "elements/segment.h"
#include "elasticity/bie_elastostatic_axi3d_mode1.h"
#include "elasticity/fullspace_iso_axisymmetry_flat_unidirectional/elastic_axi3dP0_element.h"

namespace bigwham {

// Mode-I only H-kernel for axisymmetric ring P0 element.
// DOF dimension = 1 (normal/opening displacement discontinuity only).
// Returns a 1x1 influence matrix: normal traction due to opening DD.
// NOTE: Poisson's ratio is taken as zero (as in Axi3DS0-H).
template <>
inline std::vector<double>
BieElastostaticAxi3DModeI<Segment<0>, Segment<0>, ElasticKernelType::H>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {

    double G = this->elas_.shear_modulus();

    auto src_vertices = source_elt.vertices();

    double rExt = std::sqrt(src_vertices(1, 0) * src_vertices(1, 0) +
                            src_vertices(1, 1) * src_vertices(1, 1));
    double rInt = std::sqrt(src_vertices(0, 0) * src_vertices(0, 0) +
                            src_vertices(0, 1) * src_vertices(0, 1));

    auto rec_colpts = receiver_elt.collocation_points();
    double rObs = std::sqrt(rec_colpts(0, 0) * rec_colpts(0, 0) +
                            rec_colpts(0, 1) * rec_colpts(0, 1));

    double IF = stress_disk_dislocation(rObs, rExt) - stress_disk_dislocation(rObs, rInt);

    // normal traction due to normal DD only: 2G * IF (2G = E when nu = 0)
    std::vector<double> stnl(1, 2.0 * G * IF);
    return stnl;
}

template class BieElastostaticAxi3DModeI<Segment<0>, Segment<0>, ElasticKernelType::H>;

} // namespace bigwham
