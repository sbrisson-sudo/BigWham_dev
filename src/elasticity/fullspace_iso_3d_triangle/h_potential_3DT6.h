//
// This file part of BigWham
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef BIGWHAM_H_POTENTIAL_3DT6_H
#define BIGWHAM_H_POTENTIAL_3DT6_H

#include <complex>

#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>

namespace bigwham {

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> s_ij_gen_h(double nu, std::complex<double> eix,double h, std::complex<double> d);

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> s_ij_red_h(double nu, std::complex<double> eix,double h);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_h(double nu, std::complex<double> eix,std::complex<double> d);

}
#endif //BIGWHAM_H_POTENTIAL_3DT6_H
