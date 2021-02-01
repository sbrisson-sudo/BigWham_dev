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
#include <src/core/FaceData.h>
#include <src/core/ElasticProperties.h>


namespace bie{
il::StaticArray2D<double, 3, 6> StressesKernelR0(
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) ;

bool is_stress_singular_at_given_location(double&x, double& y, double& z, double& a, double& b, bool verbose = true) ;


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

il::Array<double> point_stress_3DR0(
            il::Array<double> &observ_pt, // coordinates
            FaceData &elem_data_s, // source element
            il::Array<double> &dd, // dislacement discontinuities components
            ElasticProperties const &elas_); // elastic properties

il::Array<double> point_displacement_3DR0(
            il::Array<double> &observ_pt, // coordinates
            FaceData &elem_data_s, // source element
            il::Array<double> &dd, // dislacement discontinuities components
            ElasticProperties const &elas_); // elastic properties

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

// special cases:

double Ip33_lim_z_to_0_and_x_to_a(double& x, double& y, double& a, double& b);

double Ip33_lim_z_to_0_and_y_to_b(double& x, double& y, double& a, double& b);

}
#endif //HFPX3D_R0_ELEMENT_H