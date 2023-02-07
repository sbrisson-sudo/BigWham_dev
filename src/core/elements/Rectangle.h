//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 27.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_RECTANGLE_H
#define BIGWHAM_RECTANGLE_H
#pragma once

#include "core/BoundaryElement.h"

namespace bie{

    ///////////////////////////////////////////////////////////////////////////////////////////////
//// class for Rectangular element
    template<int p>
    class Rectangle : public BoundaryElement<3, 4, p> {
    private:
        int spatial_dimension_ = 3;

    protected:
        int n_vertices_ = 4;
        int n_nodes_ = 1;
        il::StaticArray2D<double, 4, 3> vertices_;
        double area_;

    public:
        Rectangle() ;
        ~Rectangle();

        void setElement(il::Array2D<double> xv) {
            IL_ASSERT(xv.size(0)==n_vertices_);
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
            // check if this is indeed a Rectangle !
            auto dot_1_2=il::dot(this->s_,this->t_);
            IL_EXPECT_FAST(dot_1_2==0);
            double size_s = il::norm(this->s_, il::Norm::L2);
            double size_t = il::norm(this->t_, il::Norm::L2);
            area_= size_s*size_t; //area of rectangle
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
            this->setCollocationPoints();
            this->setNodes();
        }

        void setCollocationPoints(){};
        void setNodes() {
            this->nodes_ = this->collocation_points_; // by default nodes = collocation points for 0 element
        };

        il::int_t getNumberOfNodes() const { return n_nodes_; };
        il::int_t getNumberOfCollocationPoints() const { return n_nodes_; };
        double area() const {return  area_;};
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
    Rectangle<p>::Rectangle() = default;

    template<int p>
    Rectangle<p>::~Rectangle() = default;


    /////////////////////////////////////////////////////////////////
    //   Rectangle 0
    template<>
    inline Rectangle<0>::Rectangle() { n_nodes_=1; area_=0; };

    template<>
    inline void Rectangle<0>::setNodes()  {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->nodes_ = col;
    }

    template<>
    inline void Rectangle<0>::setCollocationPoints()   {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->collocation_points_ = col;
    }

}
#endif //BIGWHAM_RECTANGLE_H
