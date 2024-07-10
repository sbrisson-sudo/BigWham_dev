//
// This file part of BigWham
//
// Created by Brice Lecampion on 04.02.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef BIGWHAM_ELASTIC_3DR0_ELEMENT_H
#define BIGWHAM_ELASTIC_3DR0_ELEMENT_H

#include <il/StaticArray2D.h>

namespace bigwham{
il::StaticArray2D<double, 3, 6> StressesKernelR0(
    double& xx, double& yy, double& zz, double& a, double& b, double& G,
    double& nu) ;

il::Array2D<double> DisplacementKernelR0(
        double& xx, double& yy, double& zz, double& a, double& b,
        double& nu) ;

}

#endif
