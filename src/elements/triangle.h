//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_TRIANGLE_H
#define BIGWHAM_TRIANGLE_H

#include "elements/polygon.h"

namespace bigwham {

template <int p> class Triangle : public Polygon<p> {
private:
    // collocation points' location parameter for linear element (not implemented)
  const double beta_ =1.5 * 0.166666666666667;
    // collocation points' location parameter for Quadratic element ( implemented)
  const double beta1_ = 0.35;  // 1.5 * 0.091576213509771 related to nodes at vertices
                         // (see documentation)
    // collocation points' location parameter for Quadratic element ( implemented)
    const double beta2_ = 0.35;  // 1.5 * 0.10810301816807 related to middle-edge nodes
                         // (see documentation)
public:
  Triangle() : Polygon<p>() {
    this->num_vertices_ = 3;
    switch (p) {
    case 0:
      this->num_nodes_ = 1;
      this->num_collocation_points_ = this->num_nodes_;
      break;
    case 1:
      this->num_nodes_ = this->spatial_dimension_;
      this->num_collocation_points_ = this->num_nodes_;
      break;
    case 2:
      this->num_nodes_ = 2 * this->spatial_dimension_;
      this->num_collocation_points_ = this->num_nodes_;
      break;
    }
  this->vertices_.Resize(this->num_vertices_, 3);
  this->collocation_points_.Resize(this->num_collocation_points_, 3);
  this->nodes_.Resize(this->num_nodes_, 3);
  }
  ~Triangle() {}

  double getBeta1() const { return beta1_; };
  double getBeta2() const { return beta2_; };

  virtual void SetCollocationPoints() override;
  virtual void SetNodes() override;
};
/* -------------------------------------------------------------------------- */
// Triangle2D: triangle element in a 2D space (r,z or x,y coordinates).
// Inherits from Polygon2D<p>; mirrors Triangle<p> but without 3D geometry.
/* -------------------------------------------------------------------------- */

template <int p> class Triangle2D : public Polygon2D<p> {
public:
  Triangle2D() : Polygon2D<p>() {
    this->num_vertices_ = 3;
    this->num_nodes_ = 1;
    this->num_collocation_points_ = 1;
    this->vertices_.Resize(this->num_vertices_, 2);
    this->collocation_points_.Resize(this->num_collocation_points_, 2);
    this->nodes_.Resize(this->num_nodes_, 2);
  }
  ~Triangle2D() {}

  virtual void SetCollocationPoints() override;
  virtual void SetNodes() override;
};

template <> inline void Triangle2D<0>::SetNodes() {
  il::Array2D<double> col{1, 2, 0.};
  for (il::int_t j = 0; j < 2; j++) col(0, j) = this->centroid_[j];
  this->nodes_ = col;
}

template <> inline void Triangle2D<0>::SetCollocationPoints() {
  il::Array2D<double> col{1, 2, 0.};
  for (il::int_t j = 0; j < 2; j++) col(0, j) = this->centroid_[j];
  this->collocation_points_ = col;
}

} // namespace bigwham
#endif // BIGWHAM_TRIANGLE_H
