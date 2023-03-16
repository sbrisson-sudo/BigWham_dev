//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.02.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_POLYGON_H
#define BIGWHAM_POLYGON_H

#include <il/Array2D.h>
#include <il/StaticArray2D.h>
#include <core/BoundaryElement.h>

namespace bie {

    // class for polygon element - with fixed number of vertex ;(
    template<int p>
    class Polygon : public BoundaryElement<3,p > {

    protected:

        int spatial_dimension_ = 3;
        int n_vertices_ ;
        int n_nodes_ = p;
        il::Array2D<double> vertices_;
        double area_;

    public:
        Polygon() ;
        ~Polygon();

        void setElement(il::Array2D<double> xv) {
            IL_EXPECT_FAST(xv.size(1)==spatial_dimension_);
            n_vertices_=xv.size(0);
            il::Array2D<double> ver{n_vertices_,spatial_dimension_,0.};
            vertices_=ver;
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

        virtual void setCollocationPoints(){};

        virtual void setNodes() {
            this->nodes_ = this->collocation_points_; // by default nodes = collocation points for 0 element
        };

        int getNumberOfVertices() const { return n_vertices_; };
        il::int_t getNumberOfNodes() const { return n_nodes_; };
        il::int_t getNumberOfCollocationPoints() const { return n_nodes_; };

        il::StaticArray2D<double, 3, 3> rotationMatrix() {
            il::StaticArray2D<double, 3, 3> R;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                R(0, i) = this->s_[i];
                R(1, i) = this->t_[i];
                R(2, i) = this->n_[i];
            }
            return R;
        }

        il::StaticArray2D<double, 3, 3> rotationMatrix_T() {
            il::StaticArray2D<double, 3, 3> Rt;
            for (il::int_t i = 0; i < spatial_dimension_; i++) {
                Rt(i, 0) = this->s_[i];
                Rt(i, 1) = this->t_[i];
                Rt(i, 2) = this->n_[i];
            }
            return Rt;
        }

        il::StaticArray<double,3> to_local(il::StaticArray<double,3>  x ){
            return il::dot(this->rotationMatrix(),x);
        }

        il::StaticArray<double,3> to_global(il::StaticArray<double,3>  x ){
            return il::dot(this->rotationMatrix_T(),x);
        }

    };

    template<int p>
    Polygon<p>::Polygon() = default;

    template<int p>
    Polygon<p>::~Polygon() = default;

    // polygon 0
    template<>
    inline Polygon<0>::Polygon() { n_nodes_=1; area_=0; };

    template<>
    inline void Polygon<0>::setNodes()  {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->nodes_ = col;
    }

    template<>
    inline void Polygon<0>::setCollocationPoints()   {
        // 0 order element: collocation at centroid
        il::Array2D<double> col{1, 3, 0.};
        for (il::int_t j = 0; j < this->spatial_dimension_; j++) {
            col(0, j) = this->centroid_[j] ;
        }
        this->collocation_points_ = col;
    }


}


#endif //BIGWHAM_POLYGON_H
