//
// This file part of BigWham
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFP_PLANESTRAININFINITE_H
#define HFP_PLANESTRAININFINITE_H

#pragma once

#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

#include "core/ElasticProperties.h"
#include "core/Mesh2D.h"
#include "core/SegmentData.h"

namespace bie {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(double h, double Ep,
                                                       double x, double y);

il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t i_col,
    const ElasticProperties &Elas, double ker_options);

// by nodal effect
il::StaticArray<double, 4> stresses_kernel_dp1_dd_nodal(il::int_t local_node_i,
                                                        double h, double Ep,
                                                        double x, double y);

il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_dp1_dd_nodal(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t s_col,
    il::int_t i_col, const ElasticProperties &Elas, double ker_options);

il::StaticArray<double, 3> point_stress_s2d_dp1_dd(
    il::StaticArray<double, 2> &observ_pt, SegmentData &source_elt,
    il::StaticArray<double, 4> &nodal_dd, ElasticProperties &Elas,
    double ker_options);
}

#endif  // HFP_PLANESTRAININFINITE_H
