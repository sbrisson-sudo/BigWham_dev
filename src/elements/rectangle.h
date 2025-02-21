//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 27.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE
// file for more details.
//

#ifndef BIGWHAM_RECTANGLE_H
#define BIGWHAM_RECTANGLE_H

#include "elements/polygon.h"

namespace bigwham {

///////////////////////////////////////////////////////////////////////////////////////////////
//// class for Rectangular element
template <int p> class Rectangle : public Polygon<p> {

public:
  Rectangle() : Polygon<p>() {
    this->num_vertices_ = 4;
    switch (p) {
    case 0:
      this->num_nodes_ = 1;
      this->num_collocation_points_ = this->num_nodes_;
      break;
    }
  this->vertices_.Resize(this->num_vertices_, 3);
  this->collocation_points_.Resize(this->num_collocation_points_, 3);
  this->nodes_.Resize(this->num_nodes_, 3);
  }
  ~Rectangle() {}

  virtual void SetCollocationPoints() override;
  virtual void SetNodes() override;
};


template <> inline void Rectangle<0>::SetNodes() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->nodes_ = col;
}

template <> inline void Rectangle<0>::SetCollocationPoints() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->collocation_points_ = col;
}

} // namespace bie
#endif // BIGWHAM_RECTANGLE_H
