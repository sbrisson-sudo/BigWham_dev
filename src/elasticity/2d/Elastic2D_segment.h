//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_ELASTIC2D_SEGMENT_H
#define BIGWHAM_ELASTIC2D_SEGMENT_H

namespace bie{

    il::StaticArray2D<double,2,2 > Ue_segment_0(double h, double G, double nu,double x1_o, double y1_o);

    il::StaticArray2D<double,2,3> Se_segment_0(double h, double G, double nu,double x1_o, double y1_o);

    il::StaticArray2D<double, 2, 3> We_segment_0(double h, double G, double nu, double x1_o, double y1_o);

}
#endif //BIGWHAM_ELASTIC2D_SEGMENT_H
