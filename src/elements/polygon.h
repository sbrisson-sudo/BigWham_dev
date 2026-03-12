//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_POLYGON_H
#define BIGWHAM_POLYGON_H

#include <il/linearAlgebra/dense/norm.h>

#include "elements/boundary_element.h"

namespace bigwham {

// class for polygon element - with fixed number of vertex ;(
template <int p> class Polygon : public BoundaryElement {

public:
  Polygon() : BoundaryElement(3, p) {}
  ~Polygon() {}

  virtual void SetElement(const il::Array2D<double> &coods_vertices) override;
  virtual void SetRotationMatrices() override;
  virtual void SetCollocationPoints() = 0;
  virtual void SetNodes() = 0;
  bool isPointOnBoundary(const std::array<double, 2>& xy_obs) const;
  bool isPointInPolygon(const std::array<double, 2>& xy_obs) const;
  double getTol() const { return tol_; };

private:
  double tol_;
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

  // Check and correct element orientation to ensure counter-clockwise ordering
  // Using the shoelace formula for signed area
  double signed_area = 0.0;
  for (il::int_t i = 0; i < num_vertices_; i++) {
    il::int_t next_i = (i + 1) % num_vertices_;
    signed_area += (this->vertices_(i, 0) * this->vertices_(next_i, 1) -
                    this->vertices_(next_i, 0) * this->vertices_(i, 1));
  }

  // If clockwise (negative signed area), reverse vertex order to make counter-clockwise
  if (signed_area < 0.0) {
    // Reverse the vertex array
    for (il::int_t i = 0; i < num_vertices_ / 2; i++) {
      for (il::int_t j = 0; j < spatial_dimension_; j++) {
        double temp = this->vertices_(i, j);
        this->vertices_(i, j) = this->vertices_(num_vertices_ - 1 - i, j);
        this->vertices_(num_vertices_ - 1 - i, j) = temp;
      }
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

  // Area: use shoelace formula for all polygons (already computed above)
  this->size_ = std::abs(signed_area) / 2.0;

  // normal s and t
  for (il::int_t k = 0; k < spatial_dimension_; ++k) {
    this->tangent1_[k] = this->tangent1_[k] / size_s;
    this->tangent2_[k] = this->tangent2_[k] / size_t;
    this->normal_[k] = this->normal_[k] / size_n;
  }

  // make tangent2 perpendicular to tangent1
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

  // We set the tolerance used for isPointOnBoundary and isPointInPolygon 
  double smallest_edge = 1e100;
  for (il::int_t i = 0; i < num_vertices_; i++) {

      double a_1_x = vertices_(i, 0);
      double a_1_y = vertices_(i, 1);
      double a_2_x = vertices_((i+1)%num_vertices_, 0);
      double a_2_y = vertices_((i+1)%num_vertices_, 1);

      double edge_length = std::sqrt( (a_2_x - a_1_x)*(a_2_x - a_1_x) + (a_2_y - a_1_y)*(a_2_y - a_1_y) );
      if (edge_length < smallest_edge) smallest_edge = edge_length;
  } 

  tol_ = smallest_edge * 1e-5;
}

template <int p>
bool Polygon<p>::isPointOnBoundary(const std::array<double, 2>& xy_obs) const
{
    auto point_to_segment_distance = [](
        const std::array<double, 2>& P,
        const std::array<double, 2>& A,
        const std::array<double, 2>& B)
    {
        double dx = B[0] - A[0];
        double dy = B[1] - A[1];

        if (dx == 0.0 && dy == 0.0) {
            // A and B are the same point
            dx = P[0] - A[0];
            dy = P[1] - A[1];
            return std::sqrt(dx * dx + dy * dy);
        }

        // Project P onto segment AB
        double t = ((P[0] - A[0]) * dx + (P[1] - A[1]) * dy) / (dx * dx + dy * dy);
        t = std::max(0.0, std::min(1.0, t));

        double proj_x = A[0] + t * dx;
        double proj_y = A[1] + t * dy;

        double dist_x = P[0] - proj_x;
        double dist_y = P[1] - proj_y;

        return std::sqrt(dist_x * dist_x + dist_y * dist_y);
    };

    for (il::int_t i = 0; i < this->num_vertices_; ++i) {
        std::array<double, 2> A = { this->vertices_(i, 0), this->vertices_(i, 1) };
        std::array<double, 2> B = { this->vertices_((i + 1) % this->num_vertices_, 0),
                                    this->vertices_((i + 1) % this->num_vertices_, 1) };

        if (point_to_segment_distance(xy_obs, A, B) <= this->tol_) {
            return true;
        }
    }

    return false;
}

template <int p>
bool Polygon<p>::isPointInPolygon(const std::array<double, 2>& xy_obs) const
{
    auto cross_sign = [](
        const std::array<double, 2>& P,
        const std::array<double, 2>& A,
        const std::array<double, 2>& B)
    {
        // Cross product AB × AP
        return (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0]);
    };

    double prev_sign = 0.0;
    bool on_edge = false;

    for (il::int_t i = 0; i < this->num_vertices_; ++i) {
        std::array<double, 2> A = {this->vertices_(i, 0), this->vertices_(i, 1)};
        std::array<double, 2> B = {this->vertices_((i + 1) % this->num_vertices_, 0),
                                   this->vertices_((i + 1) % this->num_vertices_, 1)};

        double cross = cross_sign(xy_obs, A, B);

        if (std::abs(cross) <= this->tol_) {
            // Point is close to the edge
            on_edge = true;
        }

        // Store the sign only if it's meaningfully non-zero
        if (std::abs(cross) > this->tol_) {
            if (prev_sign == 0.0) {
                prev_sign = cross;
            } else if (prev_sign * cross < 0.0) {
                // Cross product changes sign → outside
                return false;
            }
        }
    }

    // Point is inside if all crosses have same sign or point is exactly on edge
    return on_edge || (prev_sign != 0.0);
}


} // namespace bie

#endif // BIGWHAM_POLYGON_H
