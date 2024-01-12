//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.01.24.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <il/blas.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

#include "elements/rectangle.h"
#include "elasticity/bie_elastostatic.h"

#include "elastic_3dR0_mode1Cartesian_element.h"
#include "elastic_3dR0_common.h"


namespace bie {

template <>
std::vector<double>
BieElastostatic<Rectangle<0>, Rectangle<0>, ElasticKernelType::H_mode1>::influence(
    const BoundaryElement &source_elt, il::int_t i_s,
    const BoundaryElement &receiver_elt, il::int_t i_r) const {
  //  return tractions - Hypersingular elastic kernel - Rectangle 0 element
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

  dsr = source_elt.ConvertToLocal(dsr);

  // get half length of the two perpendicular edges
  il::StaticArray<double,2> a_and_b = get_a_and_b(source_elt.vertices(), 4); // a rectangle has 4 verticies
  double a = a_and_b[0], b = a_and_b[1];

  // get stress influence coefficients - in the local coordinate system of the
  // source element
  il::StaticArray2D<double, 3, 6> stress = StressesKernelR0(dsr[0], dsr[1], dsr[2], a, b, G, nu);
  // Attention!
  // It is in the reference system of the source element both in the domain and in the codomain
  // index        ->    0    1    2    3    4    5
  // DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
  // DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
  // DDz (normal) -> | sxx, syy, szz, sxy, sxz, syz  |

  // normal vector at the receiver location in the reference system of the
  // source element
  auto nr = source_elt.ConvertToLocal(receiver_elt.normal());

  // compute traction vectors at receiver element cp due to (DD1,DD2,DD3) source
  // element in the reference system of the source element
  il::Array2D<double> DDs_to_traction_local{3, 3, 0.0}; // traction vectors
  il::Array<double> traction_temp{3, 0.0};   // temporary traction vector
  il::Array2D<double> sigma_temp{3, 3, 0.0}; // temporary stress tensor
  for (int i = 0; i < 3; ++i) { // loop over the rows of Stress, i.e., over each
                                // DD component effect
    // definition of temporary stress tensor
    // i = 0 : D1
    // i = 1 : D2
    // i = 2 : D3

    sigma_temp(0, 0) = stress(i, 0); // S11
    sigma_temp(0, 1) = stress(i, 3); // S12
    sigma_temp(0, 2) = stress(i, 4); // S13
    sigma_temp(1, 0) = stress(i, 3); // S21
    sigma_temp(1, 1) = stress(i, 1); // S22
    sigma_temp(1, 2) = stress(i, 5); // S23
    sigma_temp(2, 0) = stress(i, 4); // S31
    sigma_temp(2, 1) = stress(i, 5); // S32
    sigma_temp(2, 2) = stress(i, 2); // S33

    // compute temporary traction vector
    traction_temp = il::dot(sigma_temp, nr);
    // convert to global from local src element coordinate system
    auto traction_temp_global = source_elt.ConvertToGlobal(traction_temp);
    // convert to local rec element coordinate system
    auto traction_temp_local_rec =
        receiver_elt.ConvertToLocal(traction_temp_global);

    for (int j = 0; j < 3; ++j) {
      // fill the traction vectors at receiver element cp due to (DD1,DD2,DD3)
      // source element
      DDs_to_traction_local(j, i) = traction_temp_local_rec[j];
      // | t1/D1   t1/D2  t1/D3 |
      // | t2/D1   t2/D2  t2/D3 |
      // | t3/D1   t3/D2  t3/D3 |
      // local DD & local traction
      // both in the reference system of the source element
    }
  }

  // column major format
  std::vector<double> stnl(9, 0.);
  int k = 0;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      stnl[k] = DDs_to_traction_local(i, j);
      k++;
    }
  }
  return stnl;
}

template class BieElastostatic<Rectangle<0>, Rectangle<0>, ElasticKernelType::H>;
} // namespace bie
