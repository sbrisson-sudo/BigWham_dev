//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#pragma once
#include <src/core/elements/Triangle.h>

namespace bie{

//   Triangle 0
    template<>
    Triangle<0>::Triangle() { n_nodes_=1; };

    template<>
    void Triangle<0>::setNodes()  {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] + std::numeric_limits<double>::epsilon();
        }
        this->nodes_ = col;
    }

    template<>
    void Triangle<0>::setCollocationPoints() {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] + std::numeric_limits<double>::epsilon();
        }
        this->collocation_points_ = col;
    }

// Triangle 1
    template<>
    Triangle<1>::Triangle() { n_nodes_=spatial_dimension_; };

    template<>
    void Triangle<1>::setNodes() {
        IL_EXPECT_FAST(this->spatial_dimension_ == 3 && (this->vertices_).size(0) == 3);
        il::Array2D<double> nodes{this->spatial_dimension_, this->spatial_dimension_, 0};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            for (il::int_t i = 0; j < this->spatial_dimension_; i++) {
                nodes(i, j) = this->vertices_(i, j);
            }
        }
        this->nodes_ = nodes;
    };

    template<>
    void Triangle<1>::setCollocationPoints() {
        // 1 order element: collocation points
        //   col points located on the line from centroid to vertices with offset beta from vertices
        il::Array2D<double> col{3, 3, 0.};
        // loop over collocation points
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
                col(i, j) = this->beta_ * (this->centroid_[j] - this->vertices_(i, j)) + this->vertices_(i, j);
            }
        }
        this->collocation_points_ = col;
    }

// Triangle 2
    template<>
    Triangle<2>::Triangle() { n_nodes_=2*spatial_dimension_; };

    template<>
    void Triangle<2>::setNodes() {
        // nodes for the T2
//                             0
//                           /   \
//                          /     \
//                         5       4
//                        /         \
//                       /           \
//                      1------3------2
        IL_EXPECT_FAST(this->spatial_dimension_ == 3 && (this->vertices_).size(0) == 3);
        il::Array2D<double> nodes{2 * (this->spatial_dimension_), this->spatial_dimension_, 0};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
                nodes(i, j) = this->vertices_(i, j);
            }
        }
        // 3,4,5 !! what the heck with that positioning ? why 3 is between 1 and 2 and not between 0 and 1 ?
        for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
            nodes(3, i) = (this->vertices_(1, i) + this->vertices_(2, i)) / 2.;
            nodes(4, i) = (this->vertices_(0, i) + this->vertices_(2, i)) / 2.;
            nodes(5, i) = (this->vertices_(0, i) + this->vertices_(1, i)) / 2.;
        }
        this->nodes_ = std::move(nodes);
    };


    template<>
    void Triangle<2>::setCollocationPoints() {
        // 1 order element: collocation points
        //   col points located on the line from centroid to vertices with offset beta from vertices
        IL_EXPECT_FAST(this->spatial_dimension_ == 3);
        il::Array2D<double> col{2 * this->spatial_dimension_, 3, 0.};
        // loop over collocation points
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            for (il::int_t i = 0; i < this->spatial_dimension_; i++) {
                col(i, j) = this->beta1_ * (this->centroid_[j] - this->vertices_(i, j)) + this->vertices_(i, j);
            }
        }
        this->setNodes(); // because we need the middle-edge nodes
        // loop over collocation points related to middle-edge nodes (from 3 -> 5)
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            for (il::int_t i = this->spatial_dimension_; i < 2 * this->spatial_dimension_; i++) {
                col(i, j) = this->beta2_ * (this->centroid_[j] - this->nodes_(i, j)) + this->nodes_(i, j);
            }
        }
        this->collocation_points_ = col;
    }

}