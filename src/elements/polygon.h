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

protected:
  double area_;

public:
  Polygon() : BoundaryElement(3, p) {}
  ~Polygon() {}

  virtual void set_element(const il::Array2D<double> &coods_vertices) override;
  virtual void set_rotation_matrix() override;
  virtual void set_collocation_points() = 0;
  virtual void set_nodes() = 0;
};

/* -------------------------------------------------------------------------- */
// METHOD DEFINATIONS
/* -------------------------------------------------------------------------- */

template <int p> inline void Polygon<p>::set_rotation_matrix() {
  // a_local = R a_global
  // R = |s0 s1 s2|
  //     |t0 t1 t2|
  //     |n0 n1 n2|
  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_(0, i) = this->s_[i];
    rotation_matrix_(1, i) = this->t_[i];
    rotation_matrix_(2, i) = this->n_[i];
  }

  // a_global = R.T a_local
  // R.T = |s0 t0 n0|
  //       |s1 t1 n1|
  //       |s2 t2 n2|
  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_T_(i, 0) = this->s_[i];
    rotation_matrix_T_(i, 1) = this->t_[i];
    rotation_matrix_T_(i, 2) = this->n_[i];
  }
  return;
}
/* -------------------------------------------------------------------------- */

template <> inline void Polygon<0>::set_nodes() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->nodes_ = col;
}
/* -------------------------------------------------------------------------- */

template <> inline void Polygon<0>::set_collocation_points() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->collocation_points_ = col;
}
/* -------------------------------------------------------------------------- */

template <int p>
inline void Polygon<p>::set_element(const il::Array2D<double> &xv) {
  IL_EXPECT_FAST(xv.size(1) == spatial_dimension_);
  IL_EXPECT_FAST(xv.size(0) == this->num_vertices_);
  this->vertices_.Resize(num_vertices_, spatial_dimension_);
  //
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->centroid_[j] =
        0; // always reset centroid when setting the coordinates
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
    this->s_[j] = vertices_(1, j) - vertices_(0, j);
    this->t_[j] = vertices_(num_vertices_ - 1, j) - vertices_(0, j);
  }

  double size_s = il::norm(this->s_, il::Norm::L2);
  double size_t = il::norm(this->t_, il::Norm::L2);

  // normal s and t
  for (il::int_t k = 0; k < spatial_dimension_; ++k) {
    this->s_[k] = this->s_[k] / size_s;
    this->t_[k] = this->t_[k] / size_t;
  }
  // normal;
  this->n_[0] = this->s_[1] * this->t_[2] - this->s_[2] * this->t_[1];
  this->n_[1] = this->s_[2] * this->t_[0] - this->s_[0] * this->t_[2];
  this->n_[2] = this->s_[0] * this->t_[1] - this->s_[1] * this->t_[0];
  double size_n = il::norm(this->n_, il::Norm::L2);
  for (il::int_t k = 0; k < spatial_dimension_; ++k) {
    this->n_[k] = this->n_[k] / size_n;
  }
  this->t_[0] = this->n_[1] * this->s_[2] - this->n_[2] * this->s_[1];
  this->t_[1] = this->n_[2] * this->s_[0] - this->n_[0] * this->s_[2];
  this->t_[2] = this->n_[0] * this->s_[1] - this->n_[1] * this->s_[0];
  double norm = il::norm(this->t_, il::Norm::L2);
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->t_[j] = this->t_[j] / norm;
  }
  this->set_collocation_points();
  this->set_nodes();
  this->set_rotation_matrix();
}

} // namespace bie

#endif // BIGWHAM_POLYGON_H
