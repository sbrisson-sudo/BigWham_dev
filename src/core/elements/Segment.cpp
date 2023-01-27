//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <core/elements/Segment.h>

namespace bie{

    template<>
    void Segment<0>::setCollocationPoints() {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 2, 0.};
        for (il::int_t j = 0; j < 2; j++) {
            col(0, j) = this->centroid_[j];
        }
        this->collocation_points_ = col;
    }

    template<>
    void Segment<1>::setCollocationPoints() {
        // on ref element x in [-1,1]   y=0
        // collocation @ -1/sqrt(2), 1/sqrt(2)
        il::Array2D<double> x_col{2, 2, 0.0};
        il::StaticArray<double, 2> xaux;
        x_col(0, 0) = -1. / sqrt(2.);
        x_col(1, 0) = 1. / sqrt(2.);
        auto R = this->rotationMatrix();
        for (int i = 0; i < 2; ++i) {
            xaux[0] = (size_) * x_col(i, 0) / 2.;
            xaux[1] = (size_) * x_col(i, 1) / 2.;
            xaux = il::dot(R, xaux);
            x_col(i, 0) = xaux[0] + this->centroid_[0];
            x_col(i, 1) = xaux[1] + this->centroid_[1];
        }
        this->collocation_points_ = x_col;
    }


}
