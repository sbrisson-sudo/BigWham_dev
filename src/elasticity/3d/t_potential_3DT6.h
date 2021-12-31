//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 17.12.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_T_POTENTIAL_3DT6_H
#define BIGWHAM_T_POTENTIAL_3DT6_H
#include <complex>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>

namespace bie {

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> s_ij_gen_t(double nu, std::complex<double> eix,double h, std::complex<double> d);

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> s_ij_red_t(double nu, std::complex<double> eix,double h);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_t(double nu, std::complex<double> eix,std::complex<double> d);

}
#endif //BIGWHAM_T_POTENTIAL_3DT6_H
