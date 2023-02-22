//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_0_IMPLS_H
#define BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_0_IMPLS_H
#pragma once

#include <il/linearAlgebra/dense/blas/dot.h>
#include <src/elasticity/BIE_elastostatic.h>
#include <src/core/elements/Triangle.h>
#include <src/elasticity/3d/Elastic3DT0_element.h>

namespace bie{

    template<> inline std::vector<double> BIE_elastostatic<Triangle<0>,Triangle<0>,H>::influence(Triangle<0> source_elt,il::int_t i_s,Triangle<0>  receiver_elt, il::int_t i_r) const {
    //  return tractions - Hypersingular elastic kernel - Triangle 0 element
    // source_elt : element object of the source element
    // i_s : integert for the source collocation number (0,1)
    // receiver_elt  : element object of the receiver element
    // i_r : integer for the collocation number where to compute the normal and shear stress in the receiver element
    // outputs: column-major (fortran order) vector for the displacement to tractions influence matrix

        // get constitutive parameters
        double G = elas_.getG();
        double nu = elas_.getNu();

        // get coordinates receiver cp
        il::Array2D<double> el_cp_r = receiver_elt.getCollocationPoints();

        il::StaticArray<double,3> receiver_coor{0};
        for (int i=0;i<3;i++){
            receiver_coor[i]=el_cp_r(i_r,i)+ sqrt(std::numeric_limits<double>::epsilon());
        }

        // get coordinates vertices of triangular source element
        il::StaticArray2D<double,3,3>  el_vertices_s = source_elt.getVertices(); // StaticArray

        // get stress influence coefficients - in the local coordinate system of the source element
        il::StaticArray2D<double, 3, 6> Stress= StressesKernelT0(receiver_coor,el_vertices_s,G,nu);

        // normal vector at the receiver location in the reference system of the source element
        il::StaticArray<double, 3> nr =source_elt.to_local(receiver_elt.getNormal());// il::dot(R_T_source,nr);

        // compute traction vectors at receiver element cp due to (DD1,DD2,DD3) source element
        // in the reference system of the source element
        il::StaticArray2D<double,3,3> DDs_to_traction_local{0.0}; // traction vectors
        il::StaticArray<double,3> traction_temp; // temporary traction vector
        il::StaticArray2D<double,3,3> sigma_temp{0.0}; // temporary stress tensor
        for (int i = 0; i < 3; ++i)
        { //loop over the rows of Stress, i.e., over each DD component effect
            // definition of temporary stress tensor

            sigma_temp(0,0) = Stress(i,0); // S11
            sigma_temp(0,1) = Stress(i,3); // S12
            sigma_temp(0,2) = Stress(i,4); // S13
            sigma_temp(1,0) = Stress(i,3); // S21
            sigma_temp(1,1) = Stress(i,1); // S22
            sigma_temp(1,2) = Stress(i,5); // S23
            sigma_temp(2,0) = Stress(i,4); // S31
            sigma_temp(2,1) = Stress(i,5); // S32
            sigma_temp(2,2) = Stress(i,2); // S33

            // compute temporary traction vector
            traction_temp = il::dot(sigma_temp, nr);
            for (int j = 0; j < 3; ++j) {
                // fill the traction vectors at receiver element cp due to (DD1,DD2,DD3) source element
                DDs_to_traction_local(j,i) = traction_temp[j];
                // | t1/D1   t1/D2  t1/D3 |
                // | t2/D1   t2/D2  t2/D3 |
                // | t3/D1   t3/D2  t3/D3 |
                // local DD & local traction
                // both in the reference system of the source element
            }
        }

        // local-local   only for now ....
        il::StaticArray2D<double, 3, 3>  R_T_source = source_elt.rotationMatrix_T();
        auto A_rotated = il::dot(R_T_source,DDs_to_traction_local); // back to global
        il::StaticArray2D<double, 3, 3>  R_receiver = receiver_elt.rotationMatrix();
        A_rotated = il::dot(R_receiver,A_rotated); // to receiver local
        // std vector output in column major format
        std::vector<double> stnl(9,0.);
        int k=0;
        for (int j=0;j<3;j++){
            for (int i=0;i<3;i++){
                stnl[k]=A_rotated(i,j);
                k++;
            }
        }
        return stnl;
    }

}

#endif //BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_0_IMPLS_H
