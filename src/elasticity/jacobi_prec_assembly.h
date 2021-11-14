//
// This file part of BigWham
//
// Created by carlo Peruzzo on 2019-09-01.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFP_JACOBI_PREC_ASSEMBLY_H
#define HFP_JACOBI_PREC_ASSEMBLY_H

#include <src/core/ElasticProperties.h>
#include <src/core/Mesh2D.h>
#include <src/core/SegmentData.h>

#include <src/core/Mesh3D.h>

namespace bie {


typedef il::StaticArray2D<double, 2, 4> (*vKernelCall)(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t i_col,
    const ElasticProperties &Elas, double ker_options);


il::Array<double> self_influence_elastic(
    Mesh &mesh, const bie::ElasticProperties &elas, vKernelCall KernelCall,
    double ker_options);


// todo typedef for 3D to have virtual function call  and rename the function
 il::Array<double> T2_self_influence_elastic(const il::Array<il::int_t> &permutation, bie::Mesh3D &i_meshtools,
    bie::ElasticProperties &elas,
                                           il::int_t I_want_global_DD,
                                           il::int_t I_want_global_traction);


}

#endif //HFP_JACOBI_PREC_ASSEMBLY_H
