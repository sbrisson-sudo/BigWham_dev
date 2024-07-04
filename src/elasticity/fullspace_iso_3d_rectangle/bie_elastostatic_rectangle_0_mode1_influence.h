//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.01.24.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#pragma once

#include <il/blas.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include "elements/rectangle.h"

#include "elasticity/bie_elastostatic_mode1.h"
#include "elasticity/fullspace_iso_3d_rectangle/elastic_3dR0_mode1Cartesian_element.h"
#include "elasticity/fullspace_iso_3d_rectangle/elastic_3dR0_common.h"

namespace bigwham {

template <>
std::vector<double>
BieElastostaticModeI<Rectangle<0>, Rectangle<0>, ElasticKernelType::H>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {
  //  return tractions - Hypersingular elastic kernel - Rectangle 0 element - only mode 1
  // source_elt : element object of the source element
  // i_s : integert for the source collocation number (0,1)
  // receiver_elt  : element object of the receiver element
  // i_r : integer for the collocation number where to compute the normal and
  // shear stress in the receiver element outputs: column-major (fortran order)
  // vector for the displacement to tractions influence matrix

  // get constitutive parameters
  double G = this->elas_.shear_modulus();
  double nu = this->elas_.poisson_ratio();

  // get coordinates receiver cp
  auto el_cp_r = receiver_elt.collocation_points();

  // get coordinates source cp
  auto el_cp_s = source_elt.collocation_points();

  // dsr contains the component of the distance between the source and the receiver
  il::Array<double> dsr{3};
  for (int i = 0; i < 3; ++i) { dsr[i] = el_cp_r(0,i) - el_cp_s(0,i);}

  // get half length of the two perpendicular edges
  il::StaticArray<double,2> a_and_b = get_a_and_b(source_elt.vertices(), 4); // a rectangle has 4 verticies
  double a = a_and_b[0], b = a_and_b[1];

  // get stress influence coefficient - in the local coordinate system of the
  // source element
  std::vector<double> traction(1);
  traction[0] = StressesKernelR0opening(dsr[0],dsr[1],dsr[2],a, b,G,nu);

  return traction;
}

template class BieElastostaticModeI<Rectangle<0>, Rectangle<0>, ElasticKernelType::H>;

} // namespace bie
