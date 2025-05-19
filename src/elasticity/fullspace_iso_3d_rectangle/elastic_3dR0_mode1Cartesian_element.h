//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.03.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#pragma once

#ifndef BIWGHAM_ELASTIC_3DR0_MODE1CARTESIAN_ELEMENT_H
#define BIWGHAM_ELASTIC_3DR0_MODE1CARTESIAN_ELEMENT_H
#include <il/StaticArray2D.h>

namespace bigwham{

double StressesKernelR0opening(
    double& xx, double& yy, double& zz, double& a, double& b, double& G,
    double& nu) ;
}

#endif
