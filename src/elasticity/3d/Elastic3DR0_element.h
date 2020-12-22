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

il::StaticArray2D<double, 3, 6> StressesKernelRectangularP0DD(
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) ;




#endif //HFPX3D_R0_ELEMENT_H
