//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_TRIANGLE_H
#define BIGWHAM_TRIANGLE_H

#pragma once
#include <src/core/BoundaryElement.h>

namespace bie{

    ///////////////////////////////////////////////////////////////////////////////////////////////
//// class for Triangle element
    template<int p>
    class Triangle : public BoundaryElement<3, 3, p> {
    private:
        int spatial_dimension_ = 3;
        const double beta_ = 1.5 * 0.166666666666667; // collocation points' location parameter for linear
        const double beta1_ = 0.35; // 1.5 * 0.091576213509771 related to nodes at vertices (see documentation)
        const double beta2_ = 0.35; // 1.5 * 0.10810301816807 related to middle-edge nodes (see documentation)

    protected:
        int n_vertices_ = 3;
        int n_nodes_ = 1;
        il::StaticArray2D<double, 3, 3> vertices_;
        double area_;

    public:
        Triangle() ;
        ~Triangle();

        void setTriangle(il::Array2D<double> xv) {
            //
            for (il::int_t j = 0; j < spatial_dimension_; j++) {
                this->centroid_[j] = 0; // always reset centroid when setting the coordinates
                for (il::int_t i = 0; i < n_vertices_; i++) {
                    this->vertices_(i, j) = xv(i, j);
                }
            }
            for (il::int_t j = 0; j < spatial_dimension_; j++) {
                for (il::int_t i = 0; i < n_vertices_; i++) {
                    this->centroid_[j] = this->centroid_[j] + vertices_(i, j) / n_vertices_;
                }
            }
            for (il::int_t j = 0; j < spatial_dimension_; j++) {
                this->s_[j] = vertices_(1, j) - vertices_(0, j);
                this->t_[j] = vertices_(n_vertices_ - 1, j) - vertices_(0, j);
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

            this->setCollocationPoints();
            this->setNodes();
        }

        void setCollocationPoints();

        void setNodes() {
            this->nodes_ = this->collocation_points_; // by default nodes = collocation points for 0 element
        };

        il::int_t getNumberOfNodes() const { return n_nodes_; };

        il::int_t getNumberOfCollocationPoints() const { return n_nodes_; };

        il::StaticArray2D<double, 3, 3> rotationMatrix() {
            il::StaticArray2D<double, 3, 3> R;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                R(i, 0) = this->s_[i];
                R(i, 1) = this->t_[i];
                R(i, 2) = this->n_[i];
            }
            return R;
        }

        il::StaticArray2D<double, 3, 3> rotationMatrix_T() {
            il::StaticArray2D<double, 3, 3> Rt;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                Rt(0, i) = this->s_[i];
                Rt(1, i) = this->t_[i];
                Rt(2, i) = this->n_[i];
            }
            return Rt;
        }
    };

    ////////////////////////////////////////////////////////////////////
    // templated methods implementation
    template<int p>
    Triangle<p>::Triangle() = default;

    template<int p>
    Triangle<p>::~Triangle() = default;

}
#endif //BIGWHAM_TRIANGLE_H
