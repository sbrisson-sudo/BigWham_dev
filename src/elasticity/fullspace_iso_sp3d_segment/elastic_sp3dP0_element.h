//
// This file part of BigWham
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFP_SIMPLIFIED3D_H
#define HFP_SIMPLIFIED3D_H

#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

#include "core/elastic_properties.h"
#include "core/oldies/Mesh2D.h"
#include "core/oldies/SegmentData.h"

namespace bie {


// Simplified 3D kernel - piece wise constant
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(
            double a, double b,
            double G, double nu,
            double xx, double yy);

il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_s3d_dp0_dd(
    SegmentData &source_elt,
    SegmentData &receiver_elt,
    il::int_t i_col,
    const ElasticProperties &Elas,
    double ker_options);

il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_s3d_dp0_dd_nodal(
    SegmentData &source_elt,
    SegmentData &receiver_elt,
    il::int_t s_col,
    il::int_t i_col,
    const ElasticProperties &Elas,
    double ker_options);

il::StaticArray<double, 3> point_stress_s3d_dp0_dd(
    il::StaticArray<double, 2> &observ_pt,
    SegmentData &source_elt,
    il::StaticArray<double, 4> &nodal_dd, // we use only [0] and [1]
    ElasticProperties &Elas,
    double ker_options);

}

#endif // HFP_SIMPLIFIED3D_H
