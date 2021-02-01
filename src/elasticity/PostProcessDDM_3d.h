//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 10.01.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications :: Jan. 12 2021

#ifndef HFP_POSTPROCESSDDM3D_H
#define HFP_POSTPROCESSDDM3D_H

// Inclusion from Inside Loop library
#include <il/Array2D.h>

// Inclusion from the project
#include <src/core/ElasticProperties.h>
#include <src/core/Mesh3D.h>
#include <src/core/FaceData.h>

namespace bie {
    // The typedef is used to define a new type of function pointer.
    typedef il::Array<double> (*vPPrCall3D)(
            il::Array<double> &observ_pt,
            FaceData &elem_data_s, // source element
            il::Array<double> &dd,
            ElasticProperties const &elas_); // elastic properties


    il::Array2D<double> computeStresses3D(  il::Array2D<double> &observ_pts,
                                            bie::Mesh3D &mesh,
                                            bie::ElasticProperties &elas,
                                            il::Array<double> solution,
                                            vPPrCall3D PPrCall,
                                            bool are_dd_global);

    il::Array2D<double> computeDisplacements3D(   il::Array2D<double> &observ_pts,
                                                  bie::Mesh3D &mesh,
                                                  bie::ElasticProperties &elas,
                                                  il::Array<double> solution,
                                                  vPPrCall3D PPrCall,
                                                  bool are_dd_global);

}

#endif //HFP_POSTPROCESSDDM_H
