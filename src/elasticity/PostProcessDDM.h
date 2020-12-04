//
// This file is part of HFPx2D.
//
// Created by nikolski on 2/14/2018.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFP_POSTPROCESSDDM_H
#define HFP_POSTPROCESSDDM_H

// Inclusion from Inside Loop library
#include <il/Array2D.h>

// Inclusion from the project
#include <src/core/ElasticProperties.h>
#include <src/core/Mesh2D.h>
#include <src/core/SegmentData.h>

namespace bie {

    typedef il::StaticArray<double, 3> (*vPPrCall)(
            il::StaticArray<double, 2> &observ_pt,
            SegmentData &source_elt,
            il::StaticArray<double, 4> &elt_dd, // 4 DoF for either p0 or p1 elt
            ElasticProperties &Elas,
            double ker_options);

    il::Array2D<double> computeStresses2D(il::Array2D<double> &observ_pts,
                                        bie::Mesh &mesh,
                                        bie::ElasticProperties &elas,
                                        il::Array<double> solution,
                                        vPPrCall PPrCall, double ker_options);




}

#endif //HFP_POSTPROCESSDDM_H
