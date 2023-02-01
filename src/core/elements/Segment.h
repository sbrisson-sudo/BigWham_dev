//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_SEGMENT_H
#define BIGWHAM_SEGMENT_H

#include <il/Array2D.h>
#include <src/core/BoundaryElement.h>

namespace bie{

    /////////////////////////////////////////////////////////////
    //   2D segment !
    template<int p>
    class Segment : public BoundaryElement<2, 2, p> {
    private:
        int spatial_dimension_ = 2;
    protected:
        int n_vertices_ = spatial_dimension_;
        int n_nodes_ = p + 1;
        il::StaticArray2D<double, 2, 2> vertices_;
        double size_;
    public:
        Segment();
        ~Segment();

        void setElement(il::Array2D<double> xv){
            IL_ASSERT(xv.size(0)==n_vertices_);
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
            size_ = sqrt(il::dot(this->s_, this->s_));
            // normal s and t
            for (il::int_t k = 0; k < spatial_dimension_; ++k) {
                this->s_[k] = this->s_[k] / size_;
                this->t_[k] = this->t_[k] / size_;
            }
            // normal;
            this->n_[0] = -1. * this->s_[1];
            this->n_[1] = this->s_[0];
            this->setCollocationPoints();
            this->setNodes();
        };

        il::int_t getNumberOfNodes() const { return n_nodes_; };
        il::int_t getNumberOfCollocationPoints() const { return n_nodes_; };
        double getSize() const { return size_; };
        il::StaticArray<double, 2> getTangent() const { return this->s_; };

        il::StaticArray2D<double, 2, 2> rotationMatrix() {
            il::StaticArray2D<double, 2, 2> R;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                R(i, 0) = this->s_[i];
                R(i, 1) = this->n_[i];
            }
            return R;
        }

        il::StaticArray2D<double, 2, 2> rotationMatrix_T() {
            il::StaticArray2D<double, 2, 2> R;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                R(0 ,i) = this->s_[i];
                R( 1,i) = this->n_[i];
            }
            return R;
        }

        void setCollocationPoints() ;

        void setNodes() {
            this->nodes_ = this->collocation_points_; // by default nodes = collocation points for 0 element
        }

    };

    ////////////////////////////////////////////////////////////////////
    // templated methods implementation
    template<int p>
    Segment<p>::Segment() = default;

    template<int p>
    Segment<p>::~Segment() = default;

    template<int p>
    void Segment<p>::setCollocationPoints() {}  ;

    // Zero order segment
    template<>
    inline void Segment<0>::setCollocationPoints() {
        // 0 order element: collocation at centroid
        //std::cout << " in set collo ! \n";
        il::Array2D<double> col{1, 2, 0.};
        for (il::int_t j = 0; j < 2; j++) {
            col(0, j) = this->centroid_[j];
        }
        this->collocation_points_ = col;
    }

    // Linear Segment
    template<>
    inline void Segment<1>::setCollocationPoints() {
        // on ref element x in [-1,1]   y=0
        // collocation @ -1/sqrt(2), 1/sqrt(2)
        il::Array2D<double> x_col{2, 2, 0.0};
        il::StaticArray<double, 2> xaux;
        x_col(0, 0) = -1. / sqrt(2.);
        x_col(1, 0) = 1. / sqrt(2.);
        auto R = this->rotationMatrix();
        for (int i = 0; i < 2; ++i) {
            xaux[0] = (size_) * x_col(i, 0) / 2.;
            xaux[1] = (size_) * x_col(i, 1) / 2.;
            xaux = il::dot(R, xaux);
            x_col(i, 0) = xaux[0] + this->centroid_[0];
            x_col(i, 1) = xaux[1] + this->centroid_[1];
        }
        this->collocation_points_ = x_col;
    }

}

#endif //BIGWHAM_SEGMENT_H
