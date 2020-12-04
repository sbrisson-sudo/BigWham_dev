//
// HFP project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <src/core/Mesh2D.h>
#include <src/core/SegmentData.h>

namespace bie {

////////////////////////////////////////////////////////////////////////////////
//          METHODS
////////////////////////////////////////////////////////////////////////////////

il::Array2D<double> Mesh::getCollocationPoints() {
  il::Array2D<double> Xcol{(interpolation_order_ + 1), 2, 0.},
      colPoints{connectivity_.size(0) * (interpolation_order_ + 1), 2, 0};

  il::int_t j = 0;
  for (il::int_t i = 0; i < connectivity_.size(0); i++) {
    bie::SegmentData seg = this->getElementData(i);
    Xcol = seg.CollocationPoints();
    for (il::int_t j1 = 0; j1 < seg.CollocationPoints().size(0); j1++) {
      colPoints(j, 0) = Xcol(j1, 0);
      colPoints(j, 1) = Xcol(j1, 1);
      j++;
    }
  };
  return colPoints;
};
//

////////////////////////////////////////////////////////////////////////////////
// Function returning the segment characteristic object for element ne
bie::SegmentData Mesh::getElementData(const il::int_t ne) const {
  il::StaticArray2D<double, 2, 2> Xs;

  Xs(0, 0) = this->coordinates_(this->connectivity_(ne, 0), 0);
  Xs(0, 1) = this->coordinates_(this->connectivity_(ne, 0), 1);

  Xs(1, 0) = this->coordinates_(this->connectivity_(ne, 1), 0);
  Xs(1, 1) = this->coordinates_(this->connectivity_(ne, 1), 1);

  return SegmentData(Xs, this->interpolation_order_);
};

}  // namespace bie