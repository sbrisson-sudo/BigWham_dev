//
// This file part of BigWham
//
// Created by Brice Lecampion on 23.05.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#pragma once

#ifndef BIGWHAM_ELASTIC_3DT0_ELEMENT_H
#define BIGWHAM_ELASTIC_3DT0_ELEMENT_H

#include <il/StaticArray2D.h>
#include <il/StaticArray.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include "core/elastic_properties.h"

namespace bigwham {
//

il::StaticArray2D<double, 3, 6>
StressesKernelT0(il::StaticArray<double, 3> &x,il::StaticArray2D<double, 3, 3> &xv,double G,double nu);

// Fundamental displacement kernel = displacement influence coefficients
il::StaticArray2D<double, 3, 3> DisplacementKernelT0(il::StaticArray<double, 3> &x,
                                                     il::StaticArray2D<double, 3, 3> &xv,
                                                     double nu);


} // namespace bigwham

#endif // BIGWHAM_EXE_ELASTIC3DT0_ELEMENT_H
