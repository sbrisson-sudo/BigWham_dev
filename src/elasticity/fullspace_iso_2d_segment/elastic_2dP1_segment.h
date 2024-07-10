//
// This file part of BigWham
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFP_PLANESTRAININFINITE_H
#define HFP_PLANESTRAININFINITE_H


#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

#include "core/elastic_properties.h"

namespace bigwham {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(double h, double Ep,
                                                       double x, double y);
//

// by nodal effect
il::StaticArray<double, 4> We_segment_1(il::int_t local_node_i,
                                        double h, double Ep,
                                        double x, double y);


}

#endif  // HFP_PLANESTRAININFINITE_H
