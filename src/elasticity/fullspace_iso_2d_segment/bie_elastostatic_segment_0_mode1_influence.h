//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 14.05.2026.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

#include "elements/segment.h"
#include "elasticity/bie_elastostatic_mode1.h"
#include "elasticity/fullspace_iso_2d_segment/elastic_2dP0_segment.h"

namespace bigwham {

// Mode-I only H-kernel for 2D Segment P0 element.
// DOF dimension = 1 (normal/opening displacement discontinuity only).
// Returns a 1x1 influence matrix: normal traction due to opening DD.
template <>
std::vector<double>
BieElastostaticModeI<Segment<0>, Segment<0>, ElasticKernelType::H>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {
  // return normal traction - Hypersingular elastic kernel - Segment 0 element
  // Mode I only: considers only opening (normal) displacement discontinuity
  // and returns only the normal traction.
  //
  // source_elt : source element object
  // i_s        : source collocation point index
  // receiver_elt : receiver element object
  // i_r        : receiver collocation point index
  // output: 1-element vector [St(1,1)] - normal traction due to normal DD

  // switch to the frame of the source element
  il::Array<double> xe{2, 0.0};
  auto Xmid = source_elt.centroid();
  auto r_col = receiver_elt.collocation_points();
  for (int i = 0; i < 2; ++i) {
    xe[i] = r_col(i_r, i) - Xmid[i];
  }
  auto xe_local = source_elt.ConvertToLocal(xe);

  double h = source_elt.size();

  il::StaticArray2D<double, 2, 3> stress_l =
      We_segment_0(h, this->elas_.shear_modulus(), this->elas_.poisson_ratio(),
                   xe_local[0], xe_local[1]);

  // receiver normal vector in the source local coordinate system
  auto n = source_elt.ConvertToLocal(receiver_elt.normal());

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];

  std::vector<double> stnl(1, 0.);

  // normal traction due to normal DD: St(1,1)
  stnl[0] = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) +
             n2n2 * stress_l(1, 2);

  return stnl;
}

template class BieElastostaticModeI<Segment<0>, Segment<0>, ElasticKernelType::H>;

} // namespace bigwham
