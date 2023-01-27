//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 27.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#include <src/core/elements/Rectangle.h>

namespace bie{

    //   Rectangle 0
    template<>
    Rectangle<0>::Rectangle() { n_nodes_=1; area_=0; };

    template<>
    void Rectangle<0>::setNodes()  {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->nodes_ = col;
    }

    template<>
    void Rectangle<0>::setCollocationPoints()   {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->collocation_points_ = col;
    }

}