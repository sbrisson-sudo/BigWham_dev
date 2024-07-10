//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 03.02.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// implementation of  Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, artefact-free solution.
//  Geophysical Journal International
// with modifications (notably moving to the dislocation convention of positive DD in overlap)
// we append _af for artefact free
#ifndef BIGWHAM_ELASTIC_3DT0_ELEMENT_AF_H
#define BIGWHAM_ELASTIC_3DT0_ELEMENT_AF_H

#include <il/StaticArray2D.h>
#include <il/StaticArray.h>

namespace bigwham{

    // interface for test
    il::StaticArray2D<double,3,3> AngDisDisp(double x, double y, double z, double alpha, double nu);

    // interface for test
    il::StaticArray2D<double,3,3>  TDSetupD(const il::StaticArray<double, 3> &x_obs,
                                            double alpha, double nu,
                                            const il::Array<double>& TriVertex, const il::Array<double>& SideVec);

//    il::StaticArray2D<double,3,3> Displacements_TDCS(const il::StaticArray<double, 3> &x_obs,
//                                                     const il::StaticArray2D<double, 3, 3> &P,
//                                                      double nu);

    il::StaticArray2D<double,3,3> Displacements_EFCS(const il::StaticArray<double, 3> &x_obs,
                                                     const il::StaticArray2D<double, 3, 3> &P,
                                                     double nu) ;

    il::StaticArray2D<double, 3, 3> DisplacementKernelT0_af(const il::StaticArray<double, 3> &x,
                                                            const il::StaticArray2D<double, 3, 3> &xv,
                                                            double nu);

}
#endif //BIGWHAM_ELASTIC_3DT0_ELEMENT_AF_H
