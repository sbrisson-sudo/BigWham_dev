//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_POLYGON_H
#define BIGWHAM_POLYGON_H

#include "elements/boundary_element.h"
#include <il/linearAlgebra/dense/norm.h>

namespace bie {

// class for polygon element - with fixed number of vertex ;(
template <int p> class Polygon : public BoundaryElement {

public:
  Polygon() : BoundaryElement(3, p) {}
  ~Polygon() {}

  virtual void SetElement(const il::Array2D<double> &coods_vertices) override;
  virtual void SetRotationMatrices() override;
  virtual void SetCollocationPoints() = 0;
  virtual void SetNodes() = 0;
};

/* -------------------------------------------------------------------------- */
// METHOD DEFINATIONS
/* -------------------------------------------------------------------------- */

template <int p> inline void Polygon<p>::SetRotationMatrices() {
  // a_local = R a_global
  // R = |s0 s1 s2|
  //     |t0 t1 t2|
  //     |n0 n1 n2|
  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_(0, i) = this->tangent1_[i];
    rotation_matrix_(1, i) = this->tangent2_[i];
    rotation_matrix_(2, i) = this->normal_[i];
  }

  // a_global = R.T a_local
  // R.T = |s0 t0 n0|
  //       |s1 t1 n1|
  //       |s2 t2 n2|
  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_t_(i, 0) = this->tangent1_[i];
    rotation_matrix_t_(i, 1) = this->tangent2_[i];
    rotation_matrix_t_(i, 2) = this->normal_[i];
  }
  return;
}
/* -------------------------------------------------------------------------- */

template <> inline void Polygon<0>::SetNodes() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->nodes_ = col;
}
/* -------------------------------------------------------------------------- */

template <> inline void Polygon<0>::SetCollocationPoints() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->collocation_points_ = col;
}
/* -------------------------------------------------------------------------- */

template <int p>
inline void Polygon<p>::SetElement(const il::Array2D<double> &xv) {
  IL_EXPECT_FAST(xv.size(1) == spatial_dimension_);
  IL_EXPECT_FAST(xv.size(0) == this->num_vertices_);
  this->vertices_.Resize(num_vertices_, spatial_dimension_);
  //
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->centroid_[j] =0; // always reset centroid when setting the coordinates
    for (il::int_t i = 0; i < num_vertices_; i++) {
      this->vertices_(i, j) = xv(i, j);
    }
  }
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    for (il::int_t i = 0; i < num_vertices_; i++) {
      this->centroid_[j] = this->centroid_[j] + vertices_(i, j) / num_vertices_;
    }
  }
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->tangent1_[j] = vertices_(1, j) - vertices_(0, j);
    this->tangent2_[j] = vertices_(num_vertices_ - 1, j) - vertices_(0, j);
  }

  // normal: tangent1 X tangent2
  this->normal_[0] = this->tangent1_[1] * this->tangent2_[2] -
                     this->tangent1_[2] * this->tangent2_[1];
  this->normal_[1] = this->tangent1_[2] * this->tangent2_[0] -
                     this->tangent1_[0] * this->tangent2_[2];
  this->normal_[2] = this->tangent1_[0] * this->tangent2_[1] -
                     this->tangent1_[1] * this->tangent2_[0];

  double size_s = il::norm(this->tangent1_, il::Norm::L2);
  double size_t = il::norm(this->tangent2_, il::Norm::L2);
  double size_n = il::norm(this->normal_, il::Norm::L2);

  // Area
  // assuming its a paralleogram with tangent_1 and tangent_2
  if (this->num_vertices_ == 4) {
    this->size_ = size_n;
  }

  // triangle,
  if (this->num_vertices_ == 3) {
    this->size_ = size_n / 2.;
  }

  // normal s and t
  for (il::int_t k = 0; k < spatial_dimension_; ++k) {
    this->tangent1_[k] = this->tangent1_[k] / size_s;
    this->tangent2_[k] = this->tangent2_[k] / size_t;
    this->normal_[k] = this->normal_[k] / size_n;
  }

  // make tangent2 perperdicular to tangent1
  this->tangent2_[0] = this->normal_[1] * this->tangent1_[2] -
                       this->normal_[2] * this->tangent1_[1];
  this->tangent2_[1] = this->normal_[2] * this->tangent1_[0] -
                       this->normal_[0] * this->tangent1_[2];
  this->tangent2_[2] = this->normal_[0] * this->tangent1_[1] -
                       this->normal_[1] * this->tangent1_[0];
  double norm = il::norm(this->tangent2_, il::Norm::L2);
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->tangent2_[j] = this->tangent2_[j] / norm;
  }
  this->SetRotationMatrices();
  this->SetCollocationPoints();
  this->SetNodes();
}

} // namespace bie

#endif // BIGWHAM_POLYGON_H
