//
// This file part of BigWham
//
// Created by nikolski on 2/14/2018.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


//     use open MP multithreading


//#include <elasticity/2d/ElasticS3DP0_element.h>


#include "PostProcessDDM_2d.h"

namespace bie {

il::Array2D<double> computeStresses2D(il::Array2D<double>& observ_pts,
                                    bie::Mesh2D& mesh,
                                    bie::ElasticProperties& elas,
                                    il::Array<double> solution,
                                    vPPrCall PPrCall, double ker_options) {
  // Function to get stresses (ss, xy, yy) at given points
  // due to DDs of solution + given in-situ stress
  //
  // INPUTS:
  //
  // observ_pts : observation points coordinates (x, y)
  // mesh : Mesh2D object (describing the boundary element mesh)
  // elas :: elastic properties object
  // solution :: Array containing the DD on the mesh
  // ker_options : dummy argument for kernel options
  //       here for the height of the simplified 3D elt
  //
  // OUTPUTS:
  // stresses in global (reference) coordinates (xx, xy, yy)
  //
  // CONVENTIONS:
  // positive compression / positive opening DD;
  // positive tension / positive overlap DD;

  il::Array2D<double> stress_array{observ_pts.size(0), 3, 0.0};
  il::StaticArray<double, 2> observ_pt;
  il::StaticArray<double, 4> elt_DD{0.0};
  il::StaticArray<double, 3> stress_at_pt{0.0};

  il::int_t p = mesh.interpolationOrder();

  IL_EXPECT_FAST( (mesh.numberOfElts()*(p+1)*2)==solution.size() );

  // loop on all elements
  for (il::int_t e = 0; e < mesh.numberOfElts(); ++e) {
    // get characteristic of element # e
    bie::SegmentData mysege = mesh.getElementData(e);


    // vector of DD of element e (shear, normal)
    for (il::int_t i = 0; i < (p + 1); ++i) {

      elt_DD[2 * i] = solution[(2 * (p + 1) * e) + 2*i];
      elt_DD[2 * i + 1] = solution[2 * ((p + 1) * e ) + 2*i + 1];
    }

    // loop on all observation points to compute the stresses
    for (il::int_t j = 0; j < observ_pts.size(0); ++j) {
      //   get stresses at point # j
      for (il::int_t i = 0; i < 2; ++i) {
        observ_pt[i] = observ_pts(j, i);
      }

      stress_at_pt = PPrCall(observ_pt, mysege, elt_DD, elas, ker_options);

      stress_array(j, 0) += stress_at_pt[0];
      stress_array(j, 1) += stress_at_pt[1];
      stress_array(j, 2) += stress_at_pt[2];
    }
  }

  return stress_array;

}

}  // namespace bie
