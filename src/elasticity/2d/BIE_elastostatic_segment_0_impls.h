//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 23.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_0_IMPLS_H
#define BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_0_IMPLS_H
#pragma once


#include <src/elasticity/BIE_elastostatic.h>
#include <src/core/elements/Segment.h>

#include <src/elasticity/2d/Elastic2D_segment.h>

namespace bie {

    //  U kernel - pt force displacement
    template<>
    inline std::vector<double>
    BIE_elastostatic<Segment<0>, Segment<0>, U>::influence(Segment<0> source_elt, il::int_t i_s, Segment<0> receiver_elt, il::int_t i_r) const {
        //  return  - U elastic kernel - Segment 0 element - note use of simplified3D kernel here
        // source_elt : element object of the source element
        // i_s : integert for the source collocation number (0,1)
        // receiver_elt  : element object of the receiver element
        // i_r : integer for the collocation number where to compute the normal and shear stress in the receiver element
        // outputs: column-major (fortran order) vector for the displacement to traction influence matrix

        // switch to the frame of the source element....
       // il::StaticArray2D<double, 2, 2> Rt = source_elt.rotationMatrix_T();; // transpose rotation matrix
        il::StaticArray<double, 2> xe;
        auto Xmid = source_elt.getCentroid();
        auto r_col = receiver_elt.getCollocationPoints();
        for (int i = 0; i < 2; ++i) {
            xe[i] = r_col(i_r, i) - Xmid[i];
        }
        xe = source_elt.to_local(xe);

        double h = source_elt.getSize();

        il::StaticArray2D<double, 2, 2> uik = Ue_segment_0(h, this->elas_.getG(), this->elas_.getNu(), xe[0], xe[1]);
        // displacement u_k due to tractions t_i (on the source element) -

        il::StaticArray<double, 2> uk_1; uk_1[0]=uik(0,0); uk_1[1]=uik(1,0);
        il::StaticArray<double, 2> uk_2; uk_2[0]=uik(0,1); uk_2[1]=uik(1,1);
        uk_1=source_elt.to_global(uk_1);
        uk_2=source_elt.to_global(uk_2);
        // switch back to the receiver element coor
        uk_1=receiver_elt.to_local(uk_1);
        uk_2=receiver_elt.to_local(uk_2);

        std::vector<double> stnl(4, 0.); // column major storage..
        stnl[0]=uk_1[0];stnl[1]=uk_1[1];
        stnl[2]=uk_2[0];stnl[3]=uk_2[1];

        return stnl;
    }

    //      HyperSingular Kernel   H
    template<>
    inline std::vector<double>
    BIE_elastostatic<Segment<0>, Segment<0>, H>::influence(Segment<0> source_elt, il::int_t i_s, Segment<0> receiver_elt, il::int_t i_r) const {
        //  return tractions - Hypersingular elastic kernel - Segment 0 element - note use of simplified3D kernel here
        // source_elt : element object of the source element
        // i_s : integert for the source collocation number (0,1)
        // receiver_elt  : element object of the receiver element
        // i_r : integer for the collocation number where to compute the normal and shear stress in the receiver element
        // outputs: column-major (fortran order) vector for the displacement to traction influence matrix

        // switch to the frame of the source element....
        il::StaticArray<double, 2> xe;
        auto Xmid = source_elt.getCentroid();
        auto r_col = receiver_elt.getCollocationPoints();
        for (int i = 0; i < 2; ++i) {
            xe[i] = r_col(i_r, i) - Xmid[i];
        }
        xe =source_elt.to_local(xe);

        double h = source_elt.getSize();

        il::StaticArray2D<double, 2, 3> stress_l = We_segment_0(h, this->elas_.getG(), this->elas_.getNu(), xe[0], xe[1]);

        // receiver normal and tangent vector in the source local coordinates system
        il::StaticArray<double, 2> n =source_elt.to_local(receiver_elt.getNormal());
        il::StaticArray<double, 2> s =source_elt.to_local(receiver_elt.getTangent()) ;

        double n1n1 = n[0] * n[0];
        double n2n2 = n[1] * n[1];
        double n1n2 = n[0] * n[1];
        double n1s1 = n[0] * s[0];
        double n2s2 = n[1] * s[1];
        double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

        std::vector<double> stnl(4, 0.);

// shear stress in the receiver coor sys.
//  effect of shear dd
//        St(0, 0)
        stnl[0] = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) + n2s2 * stress_l(0, 2);
//  effect of normal dd
//St(0, 1)
        stnl[2] = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) + n2s2 * stress_l(1, 2);

// normal stress in the receiver coor sys.
//  effect of shear dd
//St(1, 0)
        stnl[1] = n1n1 * stress_l(0, 0) + 2. * n1n2 * stress_l(0, 1) + n2n2 * stress_l(0, 2);
//  effect of normal dd
//St(1, 1)
        stnl[3] = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) + n2n2 * stress_l(1, 2);

        return stnl;
    }






//////////////////////////////////////////////////////////////////////////////////////////////
// FOR NOW THIS IS THE SP3DP0 element !
    template<>
    inline std::vector<double>
    BIE_elastostatic_sp3d<Segment<0>, Segment<0>, H>::influence(Segment<0> source_elt, il::int_t i_s, Segment<0> receiver_elt, il::int_t i_r) const {
        //  return tractions - Hypersingular elastic kernel - Segment 0 element - note use of simplified3D kernel here
        // source_elt : element object of the source element
        // i_s : integert for the source collocation number (0,1)
        // receiver_elt  : element object of the receiver element
        // i_r : integer for the collocation number where to compute the normal and shear stress in the receiver element
        // outputs: column-major (fortran order) vector for the displacement to traction influence matrix

        // switch to the frame of the source element....

        il::StaticArray<double, 2> xe;
        auto Xmid = source_elt.getCentroid();
        auto r_col = receiver_elt.getCollocationPoints();
        for (int i = 0; i < 2; ++i) {
            xe[i] = r_col(i_r, i) - Xmid[i];
        }
        xe = source_elt.to_local(xe);

        double h = source_elt.getSize();
        double ker_options = this->kernel_properties_[0];

        il::StaticArray2D<double, 2, 3> stress_l = stresses_kernel_s3d_p0_dd(
                h / 2., ker_options / 2., this->elas_.getG(), this->elas_.getNu(), xe[0], xe[1]);

        il::StaticArray<double, 2> n = source_elt.to_local(receiver_elt.getNormal());
        il::StaticArray<double, 2> s = source_elt.to_local(receiver_elt.getTangent());

        double n1n1 = n[0] * n[0];
        double n2n2 = n[1] * n[1];
        double n1n2 = n[0] * n[1];
        double n1s1 = n[0] * s[0];
        double n2s2 = n[1] * s[1];

        double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

        std::vector<double> stnl(4, 0.);

        // shear stress
        //  effect of shear dd
//        St(0, 0)
        stnl[0] = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
                  n2s2 * stress_l(0, 2);
        //  effect of normal dd
        //St(0, 1)
        stnl[2] = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
                  n2s2 * stress_l(1, 2);

        // normal stress
        //  effect of shear dd
        //St(1, 0)
        stnl[1] = n1n1 * stress_l(0, 0) + 2. * n1n2 * stress_l(0, 1) +
                  n2n2 * stress_l(0, 2);
        //  effect of normal dd
        //St(1, 1)
        stnl[3] = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) +
                  n2n2 * stress_l(1, 2);

        return stnl;
    }

}
#endif //BIGWHAM_BIE_ELASTOSTATIC_SEGMENT_0_IMPLS_H
