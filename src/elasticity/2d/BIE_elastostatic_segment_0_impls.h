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


//
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
        il::StaticArray2D<double, 2, 2> R = source_elt.rotationMatrix();

        il::StaticArray2D<double, 2, 2> Rt = R; // transpose rotation matrix
        Rt(0, 1) = R(1, 0);
        Rt(1, 0) = R(0, 1);

        il::StaticArray<double, 2> xe;
        auto Xmid = source_elt.getCentroid();
        auto r_col = receiver_elt.getCollocationPoints();
        for (int i = 0; i < 2; ++i) {
            xe[i] = r_col(i_r, i) - Xmid[i];
        }
        xe = il::dot(Rt, xe);

        il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.getNormal());
        il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.getTangent());

        double h = source_elt.getSize();

        il::StaticArray2D<double, 2, 3> stress_l = We_segment_0(h, this->elas_.getG(), this->elas_.getNu(), xe[0], xe[1]);

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
        stnl[0] = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) + n2s2 * stress_l(0, 2);
//  effect of normal dd
//St(0, 1)
        stnl[2] = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) + n2s2 * stress_l(1, 2);

// normal stress
//  effect of shear dd
//St(1, 0)
        stnl[1] = n1n1 * stress_l(0, 0) + 2. * n1n2 * stress_l(0, 1) + n2n2 * stress_l(0, 2);
//  effect of normal dd
//St(1, 1)
        stnl[3] = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) + n2n2 * stress_l(1, 2);

        return stnl;
    }


//
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
        il::StaticArray2D<double, 2, 2> R = source_elt.rotationMatrix();

        il::StaticArray2D<double, 2, 2> Rt = R; // transpose rotation matrix
        Rt(0, 1) = R(1, 0);
        Rt(1, 0) = R(0, 1);

        il::StaticArray<double, 2> xe;
        auto Xmid = source_elt.getCentroid();
        auto r_col = receiver_elt.getCollocationPoints();
        for (int i = 0; i < 2; ++i) {
            xe[i] = r_col(i_r, i) - Xmid[i];
        }
        xe = il::dot(Rt, xe);

        il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.getNormal());
        il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.getTangent());

        double h = source_elt.getSize();
        double ker_options = this->kernel_properties_[0];

        il::StaticArray2D<double, 2, 3> stress_l = stresses_kernel_s3d_p0_dd(
                h / 2., ker_options / 2., this->elas_.getG(), this->elas_.getNu(), xe[0], xe[1]);

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
