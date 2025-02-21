//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 31.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE
// file for more details.
//

#pragma once

#ifndef BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_1_INFLUENCE
#define BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_1_INFLUENCE


#include "elasticity/bie_elastostatic.h"
#include "elements/boundary_element.h"
#include "elements/segment.h"
#include "elements/point.h"

#include "elastic_2dP1_segment.h"
/* -------------------------------------------------------------------------- */

namespace bigwham {

template class BieElastostatic<Segment<1>, Segment<1>, ElasticKernelType::H>;
template class BieElastostatic<Segment<1>, Point<2>, ElasticKernelType::W>;

template <>
std::vector<double>
BieElastostatic<Segment<1>, Segment<1>, ElasticKernelType::H>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {
  //  return tractions - Hypersingular elastic kernel - Segment 1 element
  // source_elt : element object of the source element
  // i_s : integert for the source collocation number (0,1)
  // receiver_elt  : element object of the receiver element
  // i_r : integer for the collocation number where to compute the normal and
  // shear stress in the receiver element outputs: column-major (fortran order)
  // vector for the displacement to traction influence matrix

  // OUTPUT
  // flatten std::vector<double> of length 4 with shear traction due to shear
  // dd, normal traction due to shear dd,
  //                normal traction due to shear dd, normal traction due to
  //                normal dd,
  // old output was:
  // 2*2 matrix containing the effect of node i_s of the source element, on
  // i_r
  // first row shear traction  (effect of shear and opening DD)
  // second row normal traction (effect of shear and opening DD)

  il::Array<double> xe{2, 0.};
  auto Xmid = source_elt.centroid();
  auto r_col = receiver_elt.collocation_points();
  for (int i = 0; i < 2; ++i) {
    xe[i] = r_col(i_r, i) - Xmid[i];
  }
  xe = source_elt.ConvertToLocal(xe);

  auto n = source_elt.ConvertToLocal(receiver_elt.normal());
  auto s = source_elt.ConvertToLocal(receiver_elt.tangent1());

  double h = source_elt.size();

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  il::StaticArray<double, 4> stress_l = We_segment_1(
          i_s, h, this->elas_.young_modulus_plane_strain(), xe[0], xe[1]);

  // shear traction
  // shear dd
  double shs =
      n1s1 * stress_l[0] + n1s2pn2s1 * stress_l[1] + n2s2 * stress_l[2];
  // normal dd
  double shn =
      n1s1 * stress_l[1] + n1s2pn2s1 * stress_l[2] + n2s2 * stress_l[3];

  // normal traction
  // shear dd
  double sns = n1n1 * stress_l[0] + 2 * n1n2 * stress_l[1] + n2n2 * stress_l[2];
  // normal dd
  double snn = n1n1 * stress_l[1] + 2 * n1n2 * stress_l[2] + n2n2 * stress_l[3];

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  std::vector<double> stnl(4, 0.);
  //        il::StaticArray2D<double, 2, 2> St;

  // shear stress (node 1 shear d, normal ; coordinates 2 shear dd , normal
  // dd)

  // St(0, 0)
  stnl[0] = shs;
  // St(0, 1)
  stnl[2] = shn;

  // normal stress (coordinates 1 shear dd, normal dd ; node 2 shear dd , normal
  // dd)
  // St(1, 0)
  stnl[1] = sns;
  // St(1, 1)
  stnl[3] = snn;
  return stnl;
}
/* -------------------------------------------------------------------------- */

// Kernel W
// stress at a point due to displacement over a source segment
template <>
std::vector<double>
BieElastostatic<Segment<1>, Point<2>, ElasticKernelType::W>::influence(
        const BoundaryElement &source_elt, il::int_t i_s,
        const BoundaryElement &receiver_elt, il::int_t i_r) const {
    //  return stress tensor - Hypersingular elastic kernel - Segment 0 element
    // source_elt : element object of the source element
    // i_s : integert for the source collocation number (0,1)
    // receiver_elt  : element object of the receiver element
    // i_r : integer for the collocation number where to compute the stress tensor
    // due to shear and normal unit displacement over the source segment
    // receiver element outputs: column-major (fortran order)

    // vector for the displacement to stress influence matrix

    // switch to the frame of the source element....
    il::Array<double> xe{2, 0.0};
    auto Xmid = source_elt.centroid();
    auto r_col = receiver_elt.collocation_points();
    for (int i = 0; i < 2; ++i) {
        xe[i] = r_col(i_r, i) - Xmid[i];
    }
    auto xe_local = source_elt.ConvertToLocal(xe);

    double h = source_elt.size();

//   sxxs, sxys, syys, syyn    (knowing that  sxxn and sxyn are respectively equal to sxys and syys )
    il::StaticArray<double, 4> stress_l = We_segment_1(i_s,h, this->elas_.young_modulus_plane_strain(),
                         xe_local[0], xe_local[1]);


    // switch back to global system
    il::Array2D<double > S_ut_l{2,2,0},S_un_l{2,2,0};
    S_ut_l(0,0)=stress_l[0];    S_ut_l(0,1)=stress_l[1];
    S_ut_l(1,0)=S_ut_l(0,1);    S_ut_l(1,1)=stress_l[2];

    S_un_l(0,0)=stress_l[1];    S_un_l(0,1)=stress_l[2];
    S_un_l(1,0)=stress_l[2];    S_un_l(1,1)=stress_l[3];

    auto R = source_elt.rotation_matrix();
    auto Rt= source_elt.rotation_matrix_t();

    auto S_ut_g = il::dot(R,il::dot(S_ut_l,Rt));
    auto S_un_g = il::dot(R,il::dot(S_un_l,Rt));

    std::vector<double> stnl(6, 0.);

    //  effect of shear displacement
    //  sigma_11
    stnl[0] = S_ut_g(0,0) ;
    // sigma_22
    stnl[1] = S_ut_g(1,1) ;
    // sigma_12
    stnl[2] = S_ut_g(0,1) ;

    //  effect of normal displacement
    //  sigma_11
    stnl[3] =   S_un_g(0,0) ;
    //sigma_22
    stnl[4] = S_un_g(1,1);
    // sigma_12
    stnl[5]= S_un_g(0,1);

    return stnl;
}
/* -------------------------------------------------------------------------- */


} // namespace bie

#endif // BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_1_INFLUENCE
