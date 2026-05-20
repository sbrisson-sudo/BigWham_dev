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
#include "elasticity/bie_elastostatic_mode1_sym.h"
#include "elasticity/fullspace_iso_2d_segment/elastic_2dP0_segment.h"

namespace bigwham {

// Mode-I only H-kernel for 2D Segment P0 element with symmetry about x=0.
// DOF dimension = 1 (normal/opening displacement discontinuity only).
// The mesh is defined on the x>0 half-plane; the influence accounts for both
// the actual element and its mirror image reflected about the y-axis (x=0).
// Returns a 1x1 influence matrix: normal traction due to opening DD.
template <>
std::vector<double>
BieElastostaticModeISym<Segment<0>, Segment<0>, ElasticKernelType::H>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {
  // return normal traction - Hypersingular elastic kernel - Segment 0 element
  // Mode I only with y-axis symmetry (x=0 plane of symmetry).
  // Accounts for both the actual source element (x>0) and its mirror (x<0).
  //
  // source_elt   : source element object (centroid must have x > 0)
  // i_s          : source collocation point index
  // receiver_elt : receiver element object
  // i_r          : receiver collocation point index
  // output: 1-element vector [St(1,1)] - normal traction due to normal DD

  double h = source_elt.size();
  double G = this->elas_.shear_modulus();
  double nu = this->elas_.poisson_ratio();

  // --- Contribution from the actual source element ---

  // Observer position relative to source centroid, in source local frame
  il::Array<double> xe{2, 0.0};
  auto Xmid = source_elt.centroid();
  auto r_col = receiver_elt.collocation_points();
  for (int i = 0; i < 2; ++i) {
    xe[i] = r_col(i_r, i) - Xmid[i];
  }
  auto xe_local = source_elt.ConvertToLocal(xe);

  il::StaticArray2D<double, 2, 3> stress_l =
      We_segment_0(h, G, nu, xe_local[0], xe_local[1]);

  // Receiver normal in the source local frame
  auto n = source_elt.ConvertToLocal(receiver_elt.normal());

  double n1 = n[0];
  double n2 = n[1];
  double n1n1 = n1 * n1;
  double n2n2 = n2 * n2;
  double n1n2 = n1 * n2;

  // Normal traction due to normal DD: St(1,1)
  double stnl_direct = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) +
                       n2n2 * stress_l(1, 2);

  // --- Contribution from the mirror source element (reflected about x=0) ---
  //
  // The mirror element has:
  //   centroid:  (-Xmid[0],  Xmid[1])
  //   tangent:   (-t1, t2)  where (t1, t2) = source tangent
  //   normal:    (-n_src1, n_src2)  where (n_src1, n_src2) = source normal
  //
  // Observer displacement relative to mirror centroid (global):
  //   dx_m = r_col_x - (-Xmid[0]) = r_col_x + Xmid[0]
  //   dy_m = r_col_y -   Xmid[1]
  //
  // Mirror local frame uses rotation R_mirror = [(-t1, t2); (-n_src1, n_src2)]
  //   xe_mirror_local[0] = (-t1)*(dx_m) + t2*(dy_m)
  //   xe_mirror_local[1] = (-n_src1)*(dx_m) + n_src2*(dy_m)
  //
  // Receiver normal in mirror local frame:
  //   n_mirror[0] = (-t1)*rn_global[0] + t2*rn_global[1]
  //   n_mirror[1] = (-n_src1)*rn_global[0] + n_src2*rn_global[1]

  auto src_tangent = source_elt.tangent1();   // (t1, t2) in global
  auto src_normal  = source_elt.normal();     // (n_src1, n_src2) in global
  auto rec_normal  = receiver_elt.normal();   // (rn1, rn2) in global

  double dx_m = r_col(i_r, 0) + Xmid[0];
  double dy_m = r_col(i_r, 1) - Xmid[1];

  // Observer in mirror local frame
  double xe_m0 = -src_tangent[0] * dx_m + src_tangent[1] * dy_m;
  double xe_m1 = -src_normal[0]  * dx_m + src_normal[1]  * dy_m;

  il::StaticArray2D<double, 2, 3> stress_m =
      We_segment_0(h, G, nu, xe_m0, xe_m1);

  // Receiver normal in mirror local frame
  double nm0 = -src_tangent[0] * rec_normal[0] + src_tangent[1] * rec_normal[1];
  double nm1 = -src_normal[0]  * rec_normal[0] + src_normal[1]  * rec_normal[1];

  double nm0nm0 = nm0 * nm0;
  double nm1nm1 = nm1 * nm1;
  double nm0nm1 = nm0 * nm1;

  double stnl_mirror = nm0nm0 * stress_m(1, 0) + 2. * nm0nm1 * stress_m(1, 1) +
                       nm1nm1 * stress_m(1, 2);

  std::vector<double> stnl(1, stnl_direct + stnl_mirror);
  return stnl;
}

template class BieElastostaticModeISym<Segment<0>, Segment<0>, ElasticKernelType::H>;

} // namespace bigwham
