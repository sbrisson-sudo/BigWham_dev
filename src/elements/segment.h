//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_SEGMENT_H
#define BIGWHAM_SEGMENT_H

#include <il/Array2D.h>

#include "elements/boundary_element.h"
#include "il/container/2d/Array2D.h"

namespace bie {

template <int p> class Segment : public BoundaryElement {

private:
  double size_;

protected:
  void set_nodes();
  void set_collocation_points();

public:
  Segment() : BoundaryElement(2, p) {
    this->num_vertices_ = 2;
    this->num_nodes_ = p + 1;
    this->num_collocation_points_ = this->num_nodes_;

    this->vertices_.Resize(num_vertices_, 2);
    this->collocation_points_.Resize(num_collocation_points_, 2);
    this->nodes_.Resize(num_nodes_, 2);
  }
  ~Segment() {}

  virtual void set_element(const il::Array2D<double> &vertices) override;
  virtual void set_rotation_matrix() override;

  double get_size() const { return size_; }
};
/* -------------------------------------------------------------------------- */

template <int p> inline void Segment<p>::set_rotation_matrix() {
  // rotation matrix from e_i to e'_i is by definition R_ij= e'_i . e_j
  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_(0, i) = this->s_[i];
    rotation_matrix_(1, i) = this->n_[i];
  }

  for (il::int_t i = 0; i < spatial_dimension_; i++) {
    rotation_matrix_T_(i, 0) = this->s_[i];
    rotation_matrix_T_(i, 1) = this->n_[i];
  }
}
/* -------------------------------------------------------------------------- */

template <int p>
inline void Segment<p>::set_element(const il::Array2D<double> &xv) {
  IL_ASSERT(xv.size(0) == num_vertices_);
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->centroid_[j] =
        0.0; // always reset centroid when setting the coordinates
    for (il::int_t i = 0; i < num_vertices_; i++) {
      this->vertices_(i, j) = xv(i, j);
    }
  }
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    for (il::int_t i = 0; i < num_vertices_; i++) {
      this->centroid_[j] =
          this->centroid_[j] + this->vertices_(i, j) / num_vertices_;
    }
  }
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    this->s_[j] = this->vertices_(1, j) - this->vertices_(0, j);
    this->t_[j] = this->vertices_(num_vertices_ - 1, j) - this->vertices_(0, j);
  }
  size_ = std::sqrt(il::dot(this->s_, this->s_));
  // unit s and t
  for (il::int_t k = 0; k < spatial_dimension_; ++k) {
    this->s_[k] = this->s_[k] / size_;
    this->t_[k] = this->t_[k] / size_;
  }
  // normal;
  this->n_[0] = -1. * this->s_[1];
  this->n_[1] = this->s_[0];

  this->set_rotation_matrix();
  this->set_collocation_points();
  this->set_nodes();
}
/* -------------------------------------------------------------------------- */

// Zero order segment
template <> inline void Segment<0>::set_collocation_points() {
  // 0 order element: collocation at centroid
  // std::cout << " in set collo ! \n";
  il::Array2D<double> col{1, 2, 0.};
  for (il::int_t j = 0; j < 2; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->collocation_points_ = col;
}
/* -------------------------------------------------------------------------- */

// Linear Segment
template <> inline void Segment<1>::set_collocation_points() {
  // on ref element x in [-1,1]   y=0
  // collocation @ -1/sqrt(2), 1/sqrt(2)
  il::Array2D<double> x_col{2, 2, 0.0};
  il::Array<double> xaux{2, 0.0};
  x_col(0, 0) = -1. / sqrt(2.);
  x_col(1, 0) = 1. / sqrt(2.);
  // auto Rt = this->rotationMatrix_T();
  for (int i = 0; i < 2; ++i) {
    xaux[0] = (size_)*x_col(i, 0) / 2.;
    xaux[1] = (size_)*x_col(i, 1) / 2.;
    xaux = this->convert_to_global(xaux); // il::dot(Rt, xaux);
    x_col(i, 0) = xaux[0] + this->centroid_[0];
    x_col(i, 1) = xaux[1] + this->centroid_[1];
  }
  // std::cout << x_col(0, 0) << " " << x_col(0, 1) << " " << x_col(1, 0) << " "
  //           << x_col(1, 1) << "\n";
  this->collocation_points_ = x_col;
}
/* -------------------------------------------------------------------------- */

template <int p> inline void Segment<p>::set_nodes() {
  this->nodes_ = this->collocation_points_; // by default nodes = collocation
                                            // points for 0 element
};

} // namespace bie

#endif // BIGWHAM_SEGMENT_H
