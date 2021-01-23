//
// This file is part of BigWham_exe.
//
// Created by Brice Lecampion on 23.05.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_EXE_ELASTIC3DT0_ELEMENT_H
#define BIGWHAM_EXE_ELASTIC3DT0_ELEMENT_H

#include <il/StaticArray2D.h>
#include <src/core/FaceData.h>
#include <src/core/ElasticProperties.h>

namespace bie{

    il::StaticArray2D<double, 3, 6> StressesKernelT0(
            il::StaticArray<double, 3> &x,
            il::StaticArray2D<double, 3, 3> &xv,
            double& G,
            double& nu) ;

    // generic integrals

    // By order of appearance in stress influence coefficients due to DD1

    double i5_Xi(
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &sinAlpha);

    double i7_Xi_Xi_Xi(
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &delta
    );

    double i7_Xi_Zeta_Zeta(
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &delta
    );

    double i7_Xi(
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &sinAlpha
    );

    double i7_Xi_Xi_Zeta(
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &delta
    );

    double i5_Zeta(
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &cosAlpha
    );

    double i7_Xi_Xi_Aux(
            double &eta,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta
    );

    double i5_Zeta_Zeta_Aux(
            il::StaticArray<double, 3> &L,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &cosAlpha
    );

    double i7_Xi_Zeta(
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta
    );

    double i5_Xi_Zeta(
            il::StaticArray<double, 3> &L,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &cosAlpha
    );

    // By order of appearance in stress influence coefficients due to DD2

    double i7_Zeta_Zeta_Zeta(
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &delta
    );

    double i7_Zeta(
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &cosAlpha
    );

    double i7_Zeta_Zeta_Aux(
            double &eta,
            il::StaticArray<double, 3> &cosAlpha,
            il::StaticArray<double, 3> &Lambda,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta
    );

    double i5_Xi_Xi_Aux(
            il::StaticArray<double, 3> &L,
            il::StaticArray<double, 3> &sinAlpha,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &delta,
            il::StaticArray<double, 3> &cosAlpha
    );

    // By order of appearance in stress influence coefficients due to DD3

    double i5_Aux(
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &delta
    );

    double i7_Aux(
            double &eta,
            il::StaticArray<double, 3> &q,
            il::StaticArray<double, 3> &d,
            il::StaticArray<double, 3> &D,
            il::StaticArray<double, 3> &delta
    );

    // Fundamental displacement kernel = displacement influence coefficients
    il::StaticArray2D<double, 3, 3> DisplacementKernelT0(
            il::StaticArray<double, 3> &x,
            il::StaticArray2D<double, 3, 3> &xv,
            double &nu);

    // generic (and auxiliary) integrals first appearing for displacement influence coefficients

    // By order of appearance in displacement influence coefficients due to DD1

    double i3_Xi(
            il::StaticArray<double, 3> &chi,
            il::StaticArray<double, 3> &sinAlpha
    );

    // By order of appearance in displacement influence coefficients due to DD2

    double i3_Zeta(
            il::StaticArray<double, 3> &chi,
            il::StaticArray<double, 3> &cosAlpha
    );

//    il::Array2D<double> NodeDDtriplet_to_CPtraction_influence(
//            bie::FaceData &elem_data_s, // source element
//            bie::FaceData &elem_data_r, // receiver element
//            bie::ElasticProperties const &elas_, // elastic properties
//            il::int_t I_want_global_DD,
//            il::int_t I_want_global_traction) ;
//
//    il::Array2D<double> DisplacementKernelT0(
//            double& x, double& y, double& z, double& a, double& b,
//            double& nu) ;
//
//    il::Array2D<double> NodeDDtriplet_to_CPdisplacement_influence(
//            bie::FaceData &elem_data_s, // source element
//            bie::FaceData &elem_data_r, // receiver element
//            bie::ElasticProperties const &elas_, // elastic properties
//            il::int_t I_want_global_DD,
//            il::int_t I_want_global_displacement) ;
}

#endif  // BIGWHAM_EXE_ELASTIC3DT0_ELEMENT_H
