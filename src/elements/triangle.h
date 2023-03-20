//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_TRIANGLE_H
#define BIGWHAM_TRIANGLE_H

#include "core/elements/Polygon.h"

namespace bie {

///////////////////////////////////////////////////////////////////////////////////////////////
//// class for Triangle element
template <int p> class Triangle : public Polygon<p> {
private:
  double beta_ =
      1.5 *
      0.166666666666667; // collocation points' location parameter for linear
  double beta1_ = 0.35;  // 1.5 * 0.091576213509771 related to nodes at vertices
                         // (see documentation)
  double beta2_ = 0.35;  // 1.5 * 0.10810301816807 related to middle-edge nodes
                         // (see documentation)

protected:
  int n_vertices_ = 3;
  int n_nodes_;
  il::StaticArray2D<double, 3, 3> vertices_;
  double area_{};

public:
  Triangle();
  ~Triangle();

  void setElement(il::Array2D<double> xv) {
    IL_ASSERT(xv.size(0) == n_vertices_);
    //
    for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
      for (il::int_t i = 0; i < n_vertices_; i++) {
        this->vertices_(i, j) = xv(i, j);
      }
    }

    for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
      this->centroid_[j] =
          0.0; // always reset centroid when setting the coordinates
      for (il::int_t i = 0; i < n_vertices_; i++) {
        this->centroid_[j] = this->centroid_[j] + vertices_(i, j) / n_vertices_;
      }
    }
    for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
      this->s_[j] = vertices_(1, j) - vertices_(0, j);
      this->t_[j] = vertices_(n_vertices_ - 1, j) - vertices_(0, j);
    }
    double size_s = il::norm(this->s_, il::Norm::L2);
    double size_t = il::norm(this->t_, il::Norm::L2);
    // normal s and t
    for (il::int_t k = 0; k < this->spatial_dimension_; ++k) {
      this->s_[k] = this->s_[k] / size_s;
      this->t_[k] = this->t_[k] / size_t;
    }
    // normal;
    this->n_[0] = this->s_[1] * this->t_[2] - this->s_[2] * this->t_[1];
    this->n_[1] = this->s_[2] * this->t_[0] - this->s_[0] * this->t_[2];
    this->n_[2] = this->s_[0] * this->t_[1] - this->s_[1] * this->t_[0];
    double size_n = il::norm(this->n_, il::Norm::L2);
    for (il::int_t k = 0; k < this->spatial_dimension_; ++k) {
      this->n_[k] = this->n_[k] / size_n;
    }
    this->t_[0] = this->n_[1] * this->s_[2] - this->n_[2] * this->s_[1];
    this->t_[1] = this->n_[2] * this->s_[0] - this->n_[0] * this->s_[2];
    this->t_[2] = this->n_[0] * this->s_[1] - this->n_[1] * this->s_[0];
    double norm = il::norm(this->t_, il::Norm::L2);
    for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
      this->t_[j] = this->t_[j] / norm;
    }
    this->setCollocationPoints();
    this->setNodes();
  }

  void setCollocationPoints() override;

  void setNodes() override {
    this->nodes_ = this->collocation_points_; // by default nodes = collocation
                                              // points for 0 element
  };

  int getNumberOfVertices() const override { return n_vertices_; };
  il::int_t getNumberOfNodes() const override { return n_nodes_; };
  il::int_t getNumberOfCollocationPoints() const override { return n_nodes_; };

  il::StaticArray2D<double, 3, 3> getVertices() const { return vertices_; };
};

////////////////////////////////////////////////////////////////////
// templated methods implementation
template <int p> Triangle<p>::Triangle() = default;

template <int p> Triangle<p>::~Triangle() = default;

//   Triangle 0
template <> inline Triangle<0>::Triangle() { n_nodes_ = 1; };

template <> inline void Triangle<0>::setNodes() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->nodes_ = col;
}

template <> inline void Triangle<0>::setCollocationPoints() {
  // 0 order element: collocation at centroid
  il::Array2D<double> col{1, 3, 0.};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    col(0, j) = this->centroid_[j];
  }
  this->collocation_points_ = col;
}

// Triangle 1
template <> inline Triangle<1>::Triangle() { n_nodes_ = spatial_dimension_; };

template <> inline void Triangle<1>::setNodes() {
  IL_EXPECT_FAST(this->spatial_dimension_ == 3 &&
                 (this->vertices_).size(0) == 3);
  il::Array2D<double> nodes{this->spatial_dimension_, this->spatial_dimension_,
                            0};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    for (il::int_t i = 0; j < this->spatial_dimension_; i++) {
      nodes(i, j) = this->vertices_(i, j);
    }
  }
  this->nodes_ = nodes;
};

template <> inline void Triangle<1>::setCollocationPoints() {
  // 1 order element: collocation points
  //   col points located on the line from centroid to vertices with offset beta
  //   from vertices
  il::Array2D<double> col{3, 3, 0.};
  // loop over collocation points
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
      col(i, j) = this->beta_ * (this->centroid_[j] - this->vertices_(i, j)) +
                  this->vertices_(i, j);
    }
  }
  this->collocation_points_ = col;
}

// Triangle 2
template <> inline Triangle<2>::Triangle() {
  n_nodes_ = 2 * spatial_dimension_;
};

template <> inline void Triangle<2>::setNodes() {
  // nodes for the T2
  //                             0
  //                           /   \
//                          /     \
//                         5       4
  //                        /         \
//                       /           \
//                      1------3------2
  IL_EXPECT_FAST(this->spatial_dimension_ == 3 &&
                 (this->vertices_).size(0) == 3);
  il::Array2D<double> nodes{2 * (this->spatial_dimension_),
                            this->spatial_dimension_, 0};
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
      nodes(i, j) = this->vertices_(i, j);
    }
  }
  // 3,4,5 !! what the heck with that positioning ? why 3 is between 1 and 2 and
  // not between 0 and 1 ?
  for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
    nodes(3, i) = (this->vertices_(1, i) + this->vertices_(2, i)) / 2.;
    nodes(4, i) = (this->vertices_(0, i) + this->vertices_(2, i)) / 2.;
    nodes(5, i) = (this->vertices_(0, i) + this->vertices_(1, i)) / 2.;
  }
  this->nodes_ = std::move(nodes);
};

template <> inline void Triangle<2>::setCollocationPoints() {
  // 1 order element: collocation points
  //   col points located on the line from centroid to vertices with offset beta
  //   from vertices
  IL_EXPECT_FAST(this->spatial_dimension_ == 3);
  il::Array2D<double> col{2 * this->spatial_dimension_, 3, 0.};
  // loop over collocation points
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
      col(i, j) = this->beta1_ * (this->centroid_[j] - this->vertices_(i, j)) +
                  this->vertices_(i, j);
    }
  }
  this->setNodes(); // because we need the middle-edge nodes
  // loop over collocation points related to middle-edge nodes (from 3 -> 5)
  for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
    for (il::int_t i = this->spatial_dimension_;
         i < 2 * this->spatial_dimension_; i++) {
      col(i, j) = this->beta2_ * (this->centroid_[j] - this->nodes_(i, j)) +
                  this->nodes_(i, j);
    }
  }
  this->collocation_points_ = col;
}

} // namespace bie
#endif // BIGWHAM_TRIANGLE_H
