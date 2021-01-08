//
// This file is part of HFPx3D.
//
// Created by Brice Lecampion on 04.02.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2019.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_R0_ELEMENT_H
#define HFPX3D_R0_ELEMENT_H

#include <il/StaticArray2D.h>
#include <src/core/ElasticProperties.h>

namespace bie{
il::StaticArray2D<double, 3, 6> StressesKernelR0(
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) ;

il::Array2D<double> NodeDDtriplet_to_CPtraction_influence(
            bie::FaceData &elem_data_s, // source element
            bie::FaceData &elem_data_r, // receiver element
            bie::ElasticProperties const &elas_, // elastic properties
            il::int_t I_want_global_DD,
            il::int_t I_want_global_traction) ;

il::Array2D<double> DisplacementKernelR0(
        double& x, double& y, double& z, double& a, double& b,
        double& nu) ;

il::Array2D<double> NodeDDtriplet_to_CPdisplacement_influence(
            bie::FaceData &elem_data_s, // source element
            bie::FaceData &elem_data_r, // receiver element
            bie::ElasticProperties const &elas_, // elastic properties
            il::int_t I_want_global_DD,
            il::int_t I_want_global_displacement) ;

// first order derivatives of I(x,y,z,xi,eta)

double ip1(double& x, double& y, double& z, double& xi, double& eta) ;

double ip2(double& x, double& y, double& z, double& xi, double& eta) ;

double ip3(double& x, double& y, double& z, double& xi, double& eta) ;

// second order derivatives of I(x,y,z,xi,eta)

double ip11(double& x, double& y, double& z, double& xi, double& eta) ;

double ip12(double& x, double& y, double& z, double& xi, double& eta) ;

double ip13(double& x, double& y, double& z, double& xi, double& eta) ;

double ip22(double& x, double& y, double& z, double& xi, double& eta) ;

double ip23(double& x, double& y, double& z, double& xi, double& eta) ;

double ip33(double& x, double& y, double& z, double& xi, double& eta) ;

//// third order derivatives of I(x,y,z,xi,eta)

double ip111(double& x, double& y, double& z, double& xi, double& eta) ;

double ip112(double& x, double& y, double& z, double& xi, double& eta) ;

double ip113(double& x, double& y, double& z, double& xi, double& eta) ;

double ip122(double& x, double& y, double& z, double& xi, double& eta) ;

double ip123(double& x, double& y, double& z, double& xi, double& eta) ;

double ip133(double& x, double& y, double& z, double& xi, double& eta) ;

double ip222(double& x, double& y, double& z, double& xi, double& eta) ;

double ip223(double& x, double& y, double& z, double& xi, double& eta) ;

double ip233(double& x, double& y, double& z, double& xi, double& eta) ;

double ip333(double& x, double& y, double& z, double& xi, double& eta) ;
}
#endif //HFPX3D_R0_ELEMENT_H