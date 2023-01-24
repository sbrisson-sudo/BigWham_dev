//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BOUNDARYELEMENT_H
#define BIGWHAM_BOUNDARYELEMENT_H

#pragma

#include <limits>

#include <il/math.h>
#include <il/linearAlgebra.h>

#include <il/StaticArray.h>
#include <il/Array2D.h>
#include <il/StaticArray2D.h>
#include <il/linearAlgebra/dense/norm.h>

namespace bie {

template<int dim,int nvert,int p>
class BoundaryElement {// base class for boundary element

protected:

    int spatial_dimension_ = dim;
    int interpolation_order_ = p;
    il::StaticArray2D<double,nvert,dim> vertices_;     // vertices' coordinates in global reference system -
    il::StaticArray<double, dim> centroid_{0.0};    // centroid of the element in global system of coordinates
    il::StaticArray<double, dim> n_{0.0};     // unit vector normal to element in global system of coordinates
    il::StaticArray<double, dim> s_{0.0};     // unit vector tangent to element in global system of coordinates,
    // in the direction from vertex 0 to vertex 1
    il::StaticArray<double, dim> t_{0.0};     // unit vector tangent to element in global system of coordinates,
    // orthogonal to s_ and n_ (un-used for 2D element)
    il::Array2D<double> collocation_points_; // collocation points' coordinates in global reference system
    il::Array2D<double> nodes_;        // nodes' coordinates in global reference system - size: number of nodes x dim

public:
    BoundaryElement();
    ~BoundaryElement();

    int getNumberOfVertices() const {return this->vertices_.size(0);} ;
    il::StaticArray<double, dim> getCentroid() const {return centroid_;};
    il::StaticArray<double, dim> getNormal() const {return  n_;};
    il::StaticArray<double, dim> getTangent_1() const {return  s_;};
    il::StaticArray<double, dim> getTangent_2() const {return  t_;};
    il::StaticArray2D<double,dim,dim> rotationMatrix() const {
        il::StaticArray2D<double,dim,dim> R_{0.};
    };

    il::Array2D<double> getCollocationPoints() const {return collocation_points_;};
    virtual il::int_t getNumberOfCollocationPoints() const {return collocation_points_.size(0);};

    il::Array2D<double> getNodes() const {return nodes_;};

    virtual il::int_t getNumberOfNodes() const {return nodes_.size(0);};

    il::int_t getSpatialDimension() const {return spatial_dimension_;};
    il::int_t getInterpolationOrder() const {return interpolation_order_;};
};

    template<int dim,int nvert, int p>
    BoundaryElement<dim,nvert, p>::BoundaryElement() = default;

    template<int dim,int nvert, int p>
    BoundaryElement<dim,nvert, p>::~BoundaryElement() = default;

    //   2D segment !
    template<int p>
    class Segment : public BoundaryElement<2,2,p>{
    private:
        int dim_=2;
    protected:
        int n_vertices_=2;
        int n_nodes_ = p+1;
        il::StaticArray2D<double,2,2> vertices_;
        double size_;

    public:
        Segment() ;
        ~Segment();

         void setSegment(il::Array2D<double> xv){
            //
            for (il::int_t j = 0; j < dim_; j++) {
                this->centroid_[j] =0; // always reset centroid when setting the coordinates
                for (il::int_t i = 0; i < n_vertices_; i++) {
                    this->vertices_(i,j)=xv(i,j);
                }
            }
             for (il::int_t j = 0; j < dim_; j++) {
                 for (il::int_t i = 0; i < n_vertices_; i++) {
                     this->centroid_[j] = this->centroid_[j] + vertices_(i,j) / n_vertices_;
                 }
             }
            for (il::int_t j = 0; j < dim_; j++) {
                this->s_[j] = vertices_(1, j) - vertices_(0, j);
                this->t_[j] = vertices_(n_vertices_-1, j) - vertices_(0, j);
            }
            size_ = sqrt(il::dot(this->s_, this->s_));
            // normal s and t
            for (il::int_t k = 0; k < dim_; ++k) {
                this->s_[k] = this->s_[k] / size_;
                this->t_[k] = this->t_[k] / size_;
            }
            // normal;
            this->n_[0] = -1. * this->s_[1];
            this->n_[1] = this->s_[0];

            this->setCollocationPoints();
            this->setNodes();
        }

        il::int_t getNumberOfNodes() const {return n_nodes_;};
        il::int_t getNumberOfCollocationPoints() const {return n_nodes_;};

        double getSize() const {return size_;};
        il::StaticArray<double,2> getTangent() const {return this->s_;};

        il::StaticArray2D<double,2,2> rotationMatrix(){
            il::StaticArray2D<double,2,2> R;
                for (il::int_t i = 0; i < dim_; i++) {
                    R(i,0)=this->s_[i];
                    R(i,1)=this->n_[i];
                }
            return R;
        }

        void setCollocationPoints(){
            IL_ASSERT(p==0 || p==1);
        };

        void setNodes() {
            this->nodes_=this->collocation_points_; // by default nodes = collocation points ....
        }

    };

    template<int p>
    Segment<p>::Segment() =default ;

    template<int p>
    Segment<p>::~Segment() =default ;

    template<> inline void Segment<0>::setCollocationPoints() {
       // 0 order element: collocation at centroid
         il::Array2D<double> col{1,2,0.};
         for (il::int_t j=0;j<2;j++){
             col(0,j)=this->centroid_[j];
         }
        this->collocation_points_ = col;
    }

    template<> inline void Segment<1>::setCollocationPoints() {
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


//    // the goal of this class is to replace the BoundaryElementData / FaceData  -> in order to uniform the API call to fundamenetal kernels integration on element routines
//    // and genearlize it so an element could eventually have only a single vertices and as such works for an observation point
//    template<il::int_t dim, typename T>
//    class BoundaryElementData {
//
//    private:
//
//        il::int_t interpolation_order_;    // interpolation order: 0, 1 or 2
//        il::StaticArray<double, dim> n_{0.0};     // unit vector normal to element in global system of coordinates
//        il::StaticArray<double, dim> s_{0.0};     // unit vector tangent to element in global system of coordinates,
//        // in the direction from vertex 0 to vertex 1
//        il::StaticArray<double, dim> t_{0.0};     // unit vector tangent to element in global system of coordinates,
//        // orthogonal to s_ and n_ (un-used for 2D element)
//        il::StaticArray<double, dim> centroid_{0.0};    // centroid of the element in global system of coordinates
//        il::Array2D<double> vertices_;     // vertices' coordinates in global reference system -
//        il::int_t n_vertices_;                    // Number of vertices
//
//        il::Array2D<double> nodes_;        // nodes' coordinates in global reference system - size: number of nodes x dim
//        il::Array2D<double> collocation_points_; // collocation points' coordinates in global reference system
//
//        double size_; // size or area of the element
//
//        // stored here for legacy
//        const double beta_ = 1.5 * 0.166666666666667; // collocation points' location parameter for linear
//        // elements (see documentation)
//        // collocation points' location parameters for quadratic elements
//        const double beta1_ = 0.35; // 1.5 * 0.091576213509771 related to nodes at vertices (see documentation)
//        const double beta2_ = 0.35; // 1.5 * 0.10810301816807 related to middle-edge nodes (see documentation)
//
//    public:
//
//        BoundaryElementData(il::Array2D<double> xv, il::int_t p) {
//            IL_ASSERT(dim == 3 || dim == 2); // only 2D or 3D spatial dimension !
//            IL_ASSERT(vertices_.size(1) == dim);
//            IL_ASSERT(p >= 0 && p <= 2);
//            vertices_ = xv;
//            interpolation_order_ = p;
//            n_vertices_ = vertices_.size(0);
//            // compute the centroid coordinates
//            for (il::int_t j = 0; j < dim; j++) {
//                for (il::int_t i = 0; i < n_vertices_; i++) {
//                    centroid_[j] = centroid_[j] + vertices_(i, j) / n_vertices_;
//                }
//            }
//            // would need to treat the specific case where there would be only 1 vertex !
//
//            // compute the normal to the element --> here dimension dependent
//            // prepare vectors for cross product
//            // vec01 goes from vertex 0 to vertex 1
//            // vec02 goes from vertex 0 to vertex NoV_ (vertex 2 in case of triangular element)
//            // for the 2D segment the 2 vectors are the same
//            il::StaticArray<double, dim> vec01, vec02;
//            for (il::int_t j = 0; j < dim; j++) {
//                vec01[j] = vertices_(1, j) - vertices_(0, j);
//                vec02[j] = vertices_(n_vertices_, j) - vertices_(0, j);
//            }
//            double vec01norm = sqrt(il::dot(vec01, vec01)), vec02norm = sqrt(il::dot(vec02, vec02));
//            size_=vec01norm; // size for 2D segment
//            // scale them
//            for (il::int_t k = 0; k < dim; ++k) {
//                vec01[k] = vec01[k] / vec01norm;
//                vec02[k] = vec02[k] / vec02norm;
//            }
//            s_ = vec01;
//            // check co-linearity of vec01 and vec02 ! -> note that we have to remove this for the compatibility with 2D element
//            // and that in any cases a 3D element can not have such a degenerate case
////        double dotvec12 = il::dot(vec01, vec02);
////        if (dotvec12 > 0.95 && dotvec12 < 1.05)
////        {
////            std::cout << "FaceData - ERROR! The 2nd point of the face and the last point are collinear with the 1st point.\n Implement a way to select a different vertex to define the plane of the element\n";
////        }
//            if (dim == 3) {// cross product to get normal vector
//                n_[0] = vec01[1] * vec02[2] - vec01[2] * vec02[1];
//                n_[1] = vec01[2] * vec02[0] - vec01[0] * vec02[2];
//                n_[2] = vec01[0] * vec02[1] - vec01[1] * vec02[0];
//                // normalize normal vector
//                double norm = il::norm(n_, il::Norm
//                L2);
//                for (il::int_t j = 0; j < dim; j++) {
//                    n_[j] = n_[j] / norm;
//                }
//            } else { //2D
//                n_[0] = -1. * s_[1];
//                n_[1] = s_[0];
//            }
//            // cross product between n_ and s_
//            if (dim == 3) {
//                t_[0] = n_[1] * s_[2] - n_[2] * s_[1];
//                t_[1] = n_[2] * s_[0] - n_[0] * s_[2];
//                t_[2] = n_[0] * s_[1] - n_[1] * s_[0];
//                double norm = il::norm(t_, il::Norm  L2);
//                for (il::int_t j = 0; j < dim; j++) {
//                    t_[j] = t_[j] / norm;
//                }
//            } else { //2D
//                t_ = s_;
//            }
//
//            // compute nodes -> location of where is the unknown....
//            this->computeNodes();
//            // compute collocation points
//            if (dim == 3) {
//                this->computeCollocationPoints_3D();
//            } else {
//
//                this->computeCollocationPoints_2D();
//            }
//
//        }
//// compute Nodes
//
//        void computeNodes() {
//            // inputs
//            //   -none; it requires the previous construction of an object of the class and then
//            //   uses its member variables implicitly
//            // output
//            //   -nodes' coordinates in the global reference system
//            //     for constant DD elements (p = 0), 1 node at the centroid
//            //     for linear DD elements (p = 1), 3 nodes at the vertices
//            //     for quadratic DD elements (p = 2), 3 nodes at the vertices plus 3 nodes at the middle-edge points
//            //
//            //   the most general case (p = 2) has the following convention (for lower order elements
//            //   the convention can be easily deducted):
//            //                 0
//            //               /   \
//            //              /     \
//            //             5       4
//            //            /         \
//            //           /           \
//            //          1------3------2
//            //   the output is delivered in the following 2D array form (for p = 2):
//            //     x0 y0 z0
//            //...
//
//            // determine the number of nodes per element based on interpolation order
//            switch (interpolation_order_) {
//                case 0: {
//                    number_nodes = 1;
//                    break;
//                }
//                default: {
//                    number_nodes = n_vertices_ * interpolation_order_;
//                    break;
//                }
//            }
//            // initialize the nodes dynamic array
//            il::Array2D<double> nodes{number_nodes, dim, 0.0};
//            nodes_ = nodes;
//            // compute nodes based on interpolation order
//            switch (interpolation_order_) {
//                case 0: {
//                    // one node at the centroid of the element
//                    for (il::int_t j = 0; j < dim; j++) {
//                        nodes_(0, j) = centroid_[j];
//                    }
//                    break;
//                }
//                case 1: {
//                    //  nodes at the vertices
//                    // loop over coordinates
//                    for (il::int_t j = 0; j < dim; j++) {
//                        for (il::int_t i = 0; i < n_vertices_; i++) {
//                            nodes_(i, j) = vertices_(i, j);
//                        }
//                    }
//                    break;
//                }
//                case 2: {
//                    // WOrk only for Triangle
//                    IL_ASSERT(dim == 3 && n_vertices_ == 3);
//                    // 3 nodes at the vertices plus 3 nodes at the middle-edge points
//                    // loop over coordinates
//                    for (il::int_t j = 0; j < dim; j++) {
//                        // nodes at vertex
//                        for (il::int_t i = 0; i < n_vertices_; i++) {
//                            nodes_(i, j) = vertices_(i, j);
//                        }
//                        // nodes at the middle-edge points
//                        nodes_(3, j) = (vertices_(1, j) + vertices_(2, j)) / 2.;
//                        nodes_(4, j) = (vertices_(0, j) + vertices_(2, j)) / 2.;
//                        nodes_(5, j) = (vertices_(0, j) + vertices_(1, j)) / 2.;
//                    }
//                }
//            }
//        }
//
////
//        void computeCollocationPoints_2D() {
//            IL_ASSERT(dim==2);
//            il::Array2D<double> x_col{interpolation_order_, dim, 0.};
//            il::StaticArray<double,2> xaux{0.};
//
//            switch (interpolation_order_) {
//                case 1: {  // linear DD
//                    x_col(0, 0) = -1. / sqrt(2.);
//                    x_col(1, 0) = 1. / sqrt(2.);
//                };
//                    break;
//                case 0: {
//                    x_col(0, 0) = 0.;
//                    x_col(0, 1) = 0.;
//                };
//                    break;
//                default:
//                    std::cout << "error\n";  //  error
//                    break;
//            };
//            il::Array2D<double> R= this->rotationMatrix();
//            for (int i = 0; i < interpolation_order_ + 1; ++i) {
//                xaux[0] = (size_) * x_col(i, 0) / 2.;
//                xaux[1] = (size_) * x_col(i, 1) / 2.;
//                xaux = il::dot(R, xaux);
//                x_col(i, 0) = xaux[0] + centroid_[0];
//                x_col(i, 1) = xaux[1] + centroid_[1];
//            }
//            collocation_points_ = x_col;
//        }
//
//        void computeCollocationPoints_3D() {
//            // inputs
//            //   -none; it requires the previous construction of an object of the class and then
//            //   uses its member variables implicitly
//            // output
//            //    - none it builts the collocation_points_ private array2D
//            //    coordinates in the global reference system
//            //     for constant DD elements (p = 0), 1 CP at the centroid
//            //     for linear DD elements (p = 1), 3 CPs at a "distance" beta1 from each vertex (see documentation)
//            //     for quadratic DD elements (p = 2), 3 CPs at a distance 'beta1' from each vertex, plus 3 CPs at a
//            //     distance 'beta2' from each middle-edge node (see documentation)
//            //   the most general case (p = 2) has the following convention (for lower order elements
//            //     the convention can be easily deducted):
//            //
//            //                 0
//            //               / + \
//            //              /     \
//            //             5 +   + 4
//            //            /   (+)   \
//            //           / +   +   + \
//            //          1------3------2
//            //
//            //     the numbers indicate the nodes ordering, the sign + indicate the collocation points
//            //     location that are ordered in the same way as the nodes, the sign (+) indicate the centroid location
//            //
//            //   the output is delivered in the following 2D array form (for p = 2):
//            //     x0 y0 z0 etc....
//
//            IL_ASSERT(dim == 3);
//            // determine the number of CPs per element based on interpolation order
//            il::int_t number_cp;
//            switch (interpolation_order_) {
//                case 0: {
//                    number_cp = 1;
//                    break;
//                }
//                default: {
//                    number_cp = n_vertices_ * interpolation_order_;
//                    break;
//                }
//            }
//            il::Array2D<double> cps{number_cp, 3, 0.0};
//            collocation_points_ = cps;
//
//            // compute collocation points based on interpolation order for 3D
//            switch (interpolation_order_) {
//                case 0: {
//                    if (n_vertices_ == 3) { // if 3DT0 kernel, then we "perturbate numerically' the collocation point
//                        // coordinates to get rid of the singularity related to the projection of the triangle edges
//                        // TODO: compute the limits and add them in the kernel implementation
//                        // loop over coordinates needed because collocation_points is Array2D and xc_ is Array1D
//                        for (il::int_t i = 0; i < dim; i++) {
//                            collocation_points_(0, i) = centroid_[i] + std::numeric_limits<double>::epsilon();
//                            //   2.22045e-16; // perturbation by machine epsilon
//                        }
//                    } else { // else is 3DR0 or any other BE of zero interpolation order, thus no perturbation
//                        // loop over coordinates needed because collocation_points is Array2D and xc_ is Array1D
//                        for (il::int_t i = 0; i < dim; i++) {
//                            collocation_points_(0, i) = centroid_[i];
//                        }
//                    }
//                    break;
//                }
//                    // What is below   works ONLY for Triangle element !
//                    // todo : extent to rectangular element
//                case 1: { // linear triangle
//                    IL_ASSERT(number_cp == 3);
//                    for (il::int_t j = 0; j < dim; j++) {
//                        for (il::int_t i = 0; i < number_cp; i++) {
//                            collocation_points_(i, j) = beta_ * (centroid_[j] - vertices_(i, j)) + vertices_(i, j);
//                        }
//                    }
//                    break;
//                }
//                case 2: {
//                    // quadratic triangle
//                    // loop over collocation points related to vertices (same as case 1)
//                    for (il::int_t j = 0; j < 3; j++) {
//                        for (il::int_t i = 0; i < 3; i++) {
//                            collocation_points_(i, j) = beta1_ * (centroid_[j] - vertices_(i, j)) + vertices_(i, j);
//                        }
//                    }
//                    // loop over collocation points related to middle-edge nodes (from 3 -> 5)
//                    for (il::int_t j = 0; j < 3; j++) {
//                        for (il::int_t i = 3; i < 6; i++) {
//                            collocation_points_(i, j) = beta2_ * (centroid_[j] - nodes_(i, j)) + nodes_(i, j);
//                        }
//                    }
//                }
//            }
//        }
//
//        il::Array2D<double> rotationMatrix(bool Transposed) {
//            // should work for 2D and 3D !
//            // inputs
//            //   -Transposed; if true the matrix will rotate any vector from the global to the local coordinate system
//            //                if false the matrix will rotate any vector from the local coordinate system to the global one
//            //                The DEFAULT is FALSE
//            // output
//            //   -rotation matrix based on the local orthogonal basis (s,t,n).
//            //    the rotation matrix has the following form (in the case of Transpose = true, i.e., from global to local):
//            //    sx sy sz
//            //    tx ty tz
//            //    nx ny nz
//
//            il::Array2D<double> rotM{dim, dim, 0.};
//            if (Transposed) {
//                // loop over spatial coordinates
//                for (il::int_t i = 0; i < dim; i++) {
//                    rotM(0, i) = s_[i]; // 1st row equal to s unit vector
//                    if (dim == 3) {
//                        rotM(1, i) = t_[i]; // 2nd row equal to t unit vector if 3D
//                    }
//                    rotM(dim - 1, i) = n_[i]; // last row equal to n unit vector
//                }
//            } else {
//                // loop over spatial coordinates
//                for (il::int_t i = 0; i < dim; i++) {
//                    rotM(i, 0) = s_[i]; // 1st column equal to s unit vector
//                    if (dim == 3) {
//                        rotM(i, 1) = t_[i]; // 2nd column equal to t unit vector if 3D
//                    }
//                    rotM(i, dim - 1) = n_[i]; // last column equal to n unit vector
//                }
//            }
//            return rotM;
//        }
//
//
//    }
//

};


#endif //BIGWHAM_BOUNDARYELEMENT_H
