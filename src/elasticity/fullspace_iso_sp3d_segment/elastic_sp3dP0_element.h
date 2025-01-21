//
// This file part of BigWham
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFP_SIMPLIFIED3D_H
#define HFP_SIMPLIFIED3D_H

#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

#include "core/elastic_properties.h"

namespace bigwham {

// Simplified 3D kernel - piece wise constant
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(
            double a, double b,
            double G, double nu,
            double xx, double yy);

}

#endif // HFP_SIMPLIFIED3D_H
