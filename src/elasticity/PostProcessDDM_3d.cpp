//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Jan. 12 2021


//     use open MP multithreading


//#include <elasticity/2d/ElasticS3DP0_element.h>

#include <iostream>
#include <il/linearAlgebra.h>
#include <src/elasticity/PostProcessDDM_3d.h>
#ifndef NUMBEROFTHREADS
#define NUMBEROFTHREADS 4
#endif

namespace bie {

il::Array2D<double> computeStresses3D(il::Array2D<double>& observ_pts,
                                    bie::Mesh3D& mesh,
                                    bie::ElasticProperties& elas,
                                    il::Array<double> solution,
                                    vPPrCall3D PPrCall,
                                    bool are_dd_global)
  {
      // Function to get stresses (sig_xx, sig_yy, sig_zz, sig_xy, sig_xz, sig_yz) at given points
      // due to DDs of solution
      //
      // INPUTS:
      //
      // observ_pts : observation points coordinates (x, y, z)
      // mesh : Mesh object (describing the boundary element mesh)
      // elas :: elastic properties object
      // solution :: Array containing the DD on the mesh
      //
      // OUTPUTS:
      // stresses in global (reference) coordinates (xx, xy, yy)
      //
      // CONVENTION:
      // positive tension / positive overlap DD;

      il::Array2D<double> stress_out{observ_pts.size(0), 6, 0.0};

      il::Array<double> stress_at_pt{6, 0.0};
      il::int_t numberCollPtsElt = mesh.numberCollPtsElt();
      il::int_t numberUnknownsElt = 3 * numberCollPtsElt;

      IL_EXPECT_FAST( (mesh.numberCollPts()*3)==solution.size() );

      #pragma omp parallel for num_threads(NUMBEROFTHREADS)
      // loop on all elements
      for (il::int_t e = 0; e < mesh.numberOfElts(); ++e)
      {
            // declarations are here because of the parallelization
            il::Array<double> elt_DD(numberUnknownsElt);
            il::Array<double> observ_pt(3);

            // get characteristic of element # e
            bie::FaceData mysege = mesh.getElementData(e);

            // get the vector of DD of element e (shear, shear, normal) from the global array
            for (il::int_t i = 0; i < numberUnknownsElt; ++i) {
              elt_DD[i] = solution[e * numberUnknownsElt + i];
            }

            // if the given DDs are global we need to get them to local
            if (are_dd_global){
                elt_DD = il::dot(mysege.rotationMatrix(0),elt_DD);
            }

            // loop on all observation points to compute the stresses
            for (il::int_t j = 0; j < observ_pts.size(0); ++j)
            {
                  //   get stresses at point # j
                  for (il::int_t i = 0; i < 3; ++i) {observ_pt[i] = observ_pts(j, i);}

                  //   get stresses at point # j due to DD over the element e
                  stress_at_pt = PPrCall(observ_pt, mysege, elt_DD, elas); //elt_DD are supposed to be local

                  stress_out(j, 0) += stress_at_pt[0]; // sxx
                  stress_out(j, 1) += stress_at_pt[1]; // syy
                  stress_out(j, 2) += stress_at_pt[2]; // szz
                  stress_out(j, 3) += stress_at_pt[3]; // sxy
                  stress_out(j, 4) += stress_at_pt[4]; // sxz
                  stress_out(j, 5) += stress_at_pt[5]; // syz
            } // loop ovr the observation points
      } // loop over the elements

      return stress_out;
  } // end computeStresses3D


    il::Array2D<double> computeDisplacements3D(il::Array2D<double>& observ_pts,
                                          bie::Mesh3D& mesh,
                                          bie::ElasticProperties& elas,
                                          il::Array<double> solution,
                                          vPPrCall3D PPrCall,
                                          bool are_dd_global)
    {
        // Function to get displacements (u_x, u_y, u_z) at given points
        // due to DDs of solution
        //
        // INPUTS:
        //
        // observ_pts : observation points coordinates (x, y, z)
        // mesh : Mesh object (describing the boundary element mesh)
        // elas :: elastic properties object
        // solution :: Array containing the DD on the mesh
        //
        // OUTPUTS:
        // displacements in global (reference) coordinates (xx, xy, yy)
        //
        // CONVENTION:
        // positive tension / positive overlap DD;

        il::Array2D<double> displacement_out{observ_pts.size(0), 3, 0.0};

        il::Array<double> displacement_at_pt{3, 0.0};
        il::int_t numberCollPtsElt = mesh.numberCollPtsElt();
        il::int_t numberUnknownsElt = 3 * numberCollPtsElt;

        IL_EXPECT_FAST( (mesh.numberCollPts()*3)==solution.size() );

        #pragma omp parallel for num_threads(NUMBEROFTHREADS)
        // loop on all elements
        for (il::int_t e = 0; e < mesh.numberOfElts(); ++e)
        {
            // declarations are here because of the parallelization
            il::Array<double> elt_DD{numberUnknownsElt};
            il::Array<double> observ_pt{3};

            // get characteristic of element # e
            bie::FaceData mysege = mesh.getElementData(e);

            // get the vector of DD of element e (shear, shear, normal) from the global array
            for (il::int_t i = 0; i < numberUnknownsElt; ++i) {
                elt_DD[i] = solution[e * numberUnknownsElt + i];
            }

            // if the given DDs are global we need to get them to local
            if (are_dd_global){
                elt_DD = il::dot(mysege.rotationMatrix(0),elt_DD);
            }

            // loop on all observation points to compute the stresses
            for (il::int_t j = 0; j < observ_pts.size(0); ++j)
            {
                //   get displacements at point # j
                for (il::int_t i = 0; i < 3; ++i) {observ_pt[i] = observ_pts(j, i);}

                //   get displacement at point # j due to DD over the element e
                displacement_at_pt = PPrCall(observ_pt, mysege, elt_DD, elas);

                displacement_out(j, 0) += displacement_at_pt[0];
                displacement_out(j, 1) += displacement_at_pt[1];
                displacement_out(j, 2) += displacement_at_pt[2];

            } // loop ovr the observation points
        } // loop over the elements

        return displacement_out;
    } // end displacements3D


}  // namespace bie
