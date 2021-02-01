//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.01.21 based on a work done by Alexis Sáez on 08.06.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <il/Array.h>
#include <il/linearAlgebra.h>
#include <il/StaticArray.h>
#include <il/Array2D.h>
#include "FaceData.h"
#include <iostream>
#include <limits>

namespace bie {

    //////////////////////////////////////////////////////////////////////////
    //        CONSTRUCTOR
    //////////////////////////////////////////////////////////////////////////

    FaceData::FaceData(il::Array2D<double> xv, il::int_t p) {

        // inputs
        //   -xv: vertices' coordinates in the global reference system, in the following 2D array form:
        //     x0 y0 z0
        //     x1 y1 z1
        //     x2 y2 z2
        //     ...
        //
        //     the vertices are counterclockwise ordered:
        //                 0
        //               /   \
        //              /     \
        //             /       \
        //            1- - - - -2
        //
        //   -p: interpolation order. It can be:
        //     0, for constant DD elements
        //     1, for linear DD elements
        //     2, for quadratic DD elements
        //
        // output
        //   it builds the class object by computing:
        //     -the normal (n) and tangents (s and t) unit vectors that define the local coordinate system
        //     -the centroid (xc) coordinates
        //     -the nodes' coordinates
        //     -the collocation points' coordinates
        //     -half length of the 1st edge of an element
        //     -half length of the last edge of an element

        this->vertices_ = xv;
        this->interpolation_order_ = p;
        this->NoV_ = xv.size(0);
        // check that the number of vertexes is at least 3
        if (NoV_ < 3 || xv.size(1) < 3) {std::cout << "FaceData - ERROR! The coordinates are non sufficient to define a face\n";}
        // check that the shape of xv is (number of vertex x 3) and not (3 x number of vertex)
        if (xv.size(1) > NoV_) {std::cout << "FaceData - ERROR! Please, traspose the array of vertex coordinates and retry! \n It wasn't of size (number of vertex x 3)\n";}

        // compute unit normal vector n_

        // prepare vectors for cross product
        // vec01 goes from vertex 0 to vertex 1
        // vec02 goes from vertex 0 to vertex NoV_ (vertex 2 in case of triangular element)
        il::StaticArray<double, 3> vec01, vec02;
        vec01[0] = xv(1, 0) - xv(0, 0);
        vec01[1] = xv(1, 1) - xv(0, 1);
        vec01[2] = xv(1, 2) - xv(0, 2);
        vec02[0] = xv(this->NoV_ - 1, 0) - xv(0, 0);
        vec02[1] = xv(this->NoV_ - 1, 1) - xv(0, 1);
        vec02[2] = xv(this->NoV_ - 1, 2) - xv(0, 2);

        // dot product to check the collinearity of 1st, 2nd and last point of the face
        //     -half length of the last edge of an element
        //     -half length of the 1st edge of an element
        double vec01norm=sqrt(il::dot(vec01, vec01)), vec02norm=sqrt(il::dot(vec02, vec02));

        for (il::int_t k = 0; k < vec01.size(); ++k) {
            vec01[k]=vec01[k]/vec01norm;
            vec02[k]=vec02[k]/vec02norm;}

        double dotvec12 = il::dot(vec01, vec02);
        if (dotvec12 > 0.95 && dotvec12 < 1.05)
        {
            std::cout << "FaceData - ERROR! The 2nd point of the face and the last point are collinear with the 1st point.\n Implement a way to select a different vertex to define the plane of the element\n";
        }

        // cross product to get normal vector
        il::StaticArray<double, 3> n;
        n[0] = vec01[1] * vec02[2] - vec01[2] * vec02[1];
        n[1] = vec01[2] * vec02[0] - vec01[0] * vec02[2];
        n[2] = vec01[0] * vec02[1] - vec01[1] * vec02[0];

        // normalize normal vector
        double norm;
        norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        n[0] = n[0] / norm;
        n[1] = n[1] / norm;
        n[2] = n[2] / norm;
        this->n_ = n;

        // compute unit tangent vector s_ that goes from vertex 0 to vertex 1
        il::StaticArray<double, 3> s;
        s[0] = xv(1, 0) - xv(0, 0);
        s[1] = xv(1, 1) - xv(0, 1);
        s[2] = xv(1, 2) - xv(0, 2);

        // normalize tangent vector s_
        norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        s[0] = s[0] / norm;
        s[1] = s[1] / norm;
        s[2] = s[2] / norm;
        this->s_ = s;

        // compute unit tangent vector t_ orthogonal to n_ and s_

        // cross product between n_ and s_
        il::StaticArray<double, 3> t;
        t[0] = n[1] * s[2] - n[2] * s[1];
        t[1] = n[2] * s[0] - n[0] * s[2];
        t[2] = n[0] * s[1] - n[1] * s[0];

        // normalize tangent vector t_
        norm = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
        t[0] = t[0] / norm;
        t[1] = t[1] / norm;
        t[2] = t[2] / norm;
        this->t_ = t;

        // compute centroid of the element
        this->xc_ = this->computeCentroid(xv); // go to function to see convention and other details
        // compute nodes' coordinates depending on the interpolation order
        this->nodes_ = this->computeNodes(); // go to function to see convention and other details

        // compute collocation points' coordinates depending on the interpolation order
        this->collocation_points_ = this->computeCollocationPoints(); // go to function to see convention and other details
    }

    //////////////////////////////////////////////////////////////////////////
    //        get-set functions  - i.e. public interfaces
    //////////////////////////////////////////////////////////////////////////

    il::Array2D<double> FaceData::getNodes() { // get nodes' coordinates
        return nodes_;
    }

    double FaceData::getNoV() { // get number of vertexes
        return NoV_;
    }

    il::Array2D<double> FaceData::getCollocationPoints() { // get collocation points' coordinates
        return collocation_points_;
    }

    il::Array2D<double> FaceData::getVertices() { // get vertices' coordinates
        return vertices_;
    }

    const double FaceData::getBeta1() {
        return beta1_;
    }

    const double FaceData::getBeta2() {
        return beta2_;
    }

    il::Array<double> FaceData::getNormal(){
        il::Array<double> n{3};
        for (int i = 0; i < 3; ++i) { n[i] = n_[i]; }
        return n;
    }
    //////////////////////////////////////////////////////////////////////////
    //   Methods
    //////////////////////////////////////////////////////////////////////////
    il::StaticArray<double, 3> FaceData::computeCentroid(il::Array2D<double> xv) {
        // inputs
        //   -none; it requires the previous construction of an object of the class and then
        //   uses its member variables implicitly
        // output
        //   -centroid's coordinates in the global reference system
        //
        // In the case of:
        //  - a triangular element, the centroid is the mean value of the vertex's coordinates
        //  - a quadrilateral element, the centroid is obtained from the intersection of two lines drown from 4 points.
        //    These 4 points are the centroids of all the triangles that can be generated in a quadrilateral by connecting its vertices

        il::StaticArray<double, 3> xc;
        switch (this->NoV_) {
            case 3: {
                xc[0] = (xv(0, 0) + xv(1, 0) + xv(2, 0)) / 3.;
                xc[1] = (xv(0, 1) + xv(1, 1) + xv(2, 1)) / 3.;
                xc[2] = (xv(0, 2) + xv(1, 2) + xv(2, 2)) / 3.;
                break;
            }
            case 4: {
                il::StaticArray<double, 3> x1,x2,x3,x4;
                x2[0]=(xv(0, 0) + xv(1, 0) + xv(2, 0))/ 3.;
                x2[1]=(xv(0, 1) + xv(1, 1) + xv(2, 1))/ 3.;
                x2[2]=(xv(0, 2) + xv(1, 2) + xv(2, 2))/ 3.;

                x4[0]=(xv(0, 0) + xv(2, 0) + xv(3, 0))/ 3.;
                x4[1]=(xv(0, 1) + xv(2, 1) + xv(3, 1))/ 3.;
                x4[2]=(xv(0, 2) + xv(2, 2) + xv(3, 2))/ 3.;

                x1[0]=(xv(0, 0) + xv(1, 0) + xv(3, 0))/ 3.;
                x1[1]=(xv(0, 1) + xv(1, 1) + xv(3, 1))/ 3.;
                x1[2]=(xv(0, 2) + xv(1, 2) + xv(3, 2))/ 3.;

                x3[0]=(xv(2, 0) + xv(1, 0) + xv(3, 0))/ 3.;
                x3[1]=(xv(2, 1) + xv(1, 1) + xv(3, 1))/ 3.;
                x3[2]=(xv(2, 2) + xv(1, 2) + xv(3, 2))/ 3.;

                xc=this->intersectLines(x1, x2, x3, x4);
                break;
            }
            default:{
                std::cout << "FaceData - ERROR! the computation of the centroid in the case of a number of vertexes > 4 has not been implemented!\n";
            }
        }
        return xc;
    }


    il::Array2D<double> FaceData::computeNodes() {
        // inputs
        //   -none; it requires the previous construction of an object of the class and then
        //   uses its member variables implicitly
        // output
        //   -nodes' coordinates in the global reference system
        //     for constant DD elements (p = 0), 1 node at the centroid
        //     for linear DD elements (p = 1), 3 nodes at the vertices
        //     for quadratic DD elements (p = 2), 3 nodes at the vertices plus 3 nodes at the middle-edge points
        //
        //   the most general case (p = 2) has the following convention (for lower order elements
        //   the convention can be easily deducted):
        //                 0
        //               /   \
        //              /     \
        //             5       4
        //            /         \
        //           /           \
        //          1------3------2
        //   the output is delivered in the following 2D array form (for p = 2):
        //     x0 y0 z0
        //     x1 y1 z1
        //     x2 y2 z2
        //     x3 y3 z3
        //     x4 y4 z4
        //     x5 y5 z5

        // determine the number of nodes per element based on interpolation order
        il::int_t number_nodes;
        switch (this->interpolation_order_) {
            case 0: {
                number_nodes = 1;
                break;
            }
            case 1: {
                number_nodes = 3;
                break;
            }
            case 2: {
                number_nodes = 6;
            }
        }
        // initialize the nodes dynamic array
        il::Array2D<double> nodes{number_nodes, 3, 0.0};

        // compute nodes based on interpolation order
        switch (this->interpolation_order_) {
            case 0: {
                // one node at the centroid of the element
                nodes(0, 0) = this->xc_[0];
                nodes(0, 1) = this->xc_[1];
                nodes(0, 2) = this->xc_[2];
                break;
            }
            case 1: {
                // three nodes at the vertices
                // loop over coordinates
                for (il::int_t i = 0; i < 3; i++) {
                    nodes(0, i) = this->vertices_(0, i);
                    nodes(1, i) = this->vertices_(1, i);
                    nodes(2, i) = this->vertices_(2, i);
                }
                break;
            }
            case 2: {
                // 3 nodes at the vertices plus 3 nodes at the middle-edge points
                // loop over coordinates
                for (il::int_t i = 0; i < 3; i++) {
                    // nodes at the vertices
                    nodes(0, i) = this->vertices_(0, i);
                    nodes(1, i) = this->vertices_(1, i);
                    nodes(2, i) = this->vertices_(2, i);
                    // nodes at the middle-edge points
                    nodes(3, i) = (this->vertices_(1, i) + this->vertices_(2, i)) / 2.;
                    nodes(4, i) = (this->vertices_(0, i) + this->vertices_(2, i)) / 2.;
                    nodes(5, i) = (this->vertices_(0, i) + this->vertices_(1, i)) / 2.;
                }
            }
        }
        return nodes;
    }

    il::Array2D<double> FaceData::computeCollocationPoints() {
        // inputs
        //   -none; it requires the previous construction of an object of the class and then
        //   uses its member variables implicitly
        // output
        //   -collocation points' coordinates in the global reference system
        //     for constant DD elements (p = 0), 1 CP at the centroid
        //     for linear DD elements (p = 1), 3 CPs at a "distance" beta1 from each vertex (see documentation)
        //     for quadratic DD elements (p = 2), 3 CPs at a distance 'beta1' from each vertex, plus 3 CPs at a
        //     distance 'beta2' from each middle-edge node (see documentation)

        //   the most general case (p = 2) has the following convention (for lower order elements
        //     the convention can be easily deducted):
        //
        //                 0
        //               / + \
        //              /     \
        //             5 +   + 4
        //            /   (+)   \
        //           / +   +   + \
        //          1------3------2
        //
        //     the numbers indicate the nodes ordering, the sign + indicate the collocation points
        //     location that are ordered in the same way as the nodes, the sign (+) indicate the centroid location
        //
        //   the output is delivered in the following 2D array form (for p = 2):
        //     x0 y0 z0
        //     x1 y1 z1
        //     x2 y2 z2
        //     x3 y3 z3
        //     x4 y4 z4
        //     x5 y5 z5

        // determine the number of CPs per element based on interpolation order
        il::int_t number_cp;
        switch (this->interpolation_order_) {
            case 0: {
                number_cp = 1;
                break;
            }
            case 1: {
                number_cp = 3;
                break;
            }
            case 2: {
                number_cp = 6;
            }
        }
        il::Array2D<double> cps{number_cp, 3, 0.0};

        // compute collocation points based on interpolation order

        switch (this->interpolation_order_) {
            case 0: {
                // loop over coordinates needed because collocation_points is Array2D and xc_ is Array1D
                for (il::int_t i = 0; i < 3; i++) {
                    cps(0, i) = this->xc_[i];
                }
                break;
            }
            case 1: {
                il::StaticArray<double, 3> vec_beta; // vector that goes from vertex to centroid
                // loop over collocation points
                for (il::int_t i = 0; i < 3; i++) {
                    vec_beta[0] = this->beta_ * (this->xc_[0] - this->vertices_(i, 0));
                    vec_beta[1] = this->beta_ * (this->xc_[1] - this->vertices_(i, 1));
                    vec_beta[2] = this->beta_ * (this->xc_[2] - this->vertices_(i, 2));

                    cps(i, 0) = vec_beta[0] + this->vertices_(i, 0);
                    cps(i, 1) = vec_beta[1] + this->vertices_(i, 1);
                    cps(i, 2) = vec_beta[2] + this->vertices_(i, 2);
                }
                break;
            }
            case 2: {
                il::StaticArray<double, 3> vec_beta1; // vector that goes from vertex to centroid
                // loop over collocation points related to vertices (same as case 1)
                for (il::int_t i = 0; i < 3; i++) {
                    vec_beta1[0] = this->beta1_ * (this->xc_[0] - this->vertices_(i, 0));
                    vec_beta1[1] = this->beta1_ * (this->xc_[1] - this->vertices_(i, 1));
                    vec_beta1[2] = this->beta1_ * (this->xc_[2] - this->vertices_(i, 2));

                    cps(i, 0) = vec_beta1[0] + this->vertices_(i, 0);
                    cps(i, 1) = vec_beta1[1] + this->vertices_(i, 1);
                    cps(i, 2) = vec_beta1[2] + this->vertices_(i, 2);
                }
                il::StaticArray<double, 3> vec_beta2; // vector that goes from middle-edge node to centroid
                il::Array2D<double> nodes = this->computeNodes(); // because we need the middle-edge nodes
                // loop over collocation points related to middle-edge nodes (from 3 -> 5)
                for (il::int_t i = 3; i < 6; i++) {
                    vec_beta2[0] = this->beta2_ * (this->xc_[0] - nodes(i, 0));
                    vec_beta2[1] = this->beta2_ * (this->xc_[1] - nodes(i, 1));
                    vec_beta2[2] = this->beta2_ * (this->xc_[2] - nodes(i, 2));

                    cps(i, 0) = vec_beta2[0] + nodes(i, 0);
                    cps(i, 1) = vec_beta2[1] + nodes(i, 1);
                    cps(i, 2) = vec_beta2[2] + nodes(i, 2);
                }
            }
        }
        return cps;
    }

    il::Array2D<double> FaceData::rotationMatrix(bool Transposed){
        // inputs
        //   -Transposed; if false the matrix will rotate any vector from the global to local coordinate system
        //                if true the matrix will rotate any vector from the local coordinate system to the global one
        //                The DEFAULT is FALSE
        // output
        //   -rotation matrix, from global to local, based on the local orthogonal basis (s,t,n).
        //    the rotation matrix has the following form:
        //    sx sy sz
        //    tx ty tz
        //    nx ny nz

        il::Array2D<double> rotM{3,3};
        if (Transposed){
            // loop over spatial coordinates
            for (il::int_t i = 0; i < 3; i++) {
                rotM(0,i) = this->s_[i]; // 1st row equal to s unit vector
                rotM(1,i) = this->t_[i]; // 2nd row equal to t unit vector
                rotM(2,i) = this->n_[i]; // 3rd row equal to n unit vector
            }
        }
        else{
            // loop over spatial coordinates
            for (il::int_t i = 0; i < 3; i++) {
                rotM(i,0) = this->s_[i]; // 1st column equal to s unit vector
                rotM(i,1) = this->t_[i]; // 2nd column equal to t unit vector
                rotM(i,2) = this->n_[i]; // 3rd column equal to n unit vector
            }
        }
        return rotM;
    }

    il::StaticArray<double, 3> FaceData::intersectLines(il::StaticArray<double, 3> x1,
                                                   il::StaticArray<double, 3> x2,
                                                   il::StaticArray<double, 3> x3,
                                                   il::StaticArray<double, 3> x4){
        // used to compute the centroid of a general quadrilateral
        // inputs
        //   -x1, x2, x3, x4; {x1,x3} and {x2,x4} defines two intersecting segments
        // output
        //   -intersection point

        il::StaticArray<double, 3> D31, D42, D12, D, intersectingPoint;

        for (il::int_t i = 0; i < 3; i++) {
            D31[i] = x3[i] - x1[i];
            D42[i] = x4[i] - x2[i];
            D12[i] = x1[i] - x2[i];

        }

        D[0] = D31[0] * D42[1] - D42[0] * D31[1];
        D[1] = D31[0] * D42[2] - D42[0] * D31[2];
        D[2] = D31[1] * D42[2] - D42[1] * D31[2];
        il::int_t index, i, j, k;
        index = this->maxAbsArray(D);
        switch (index) {
            case 0: { i=0; j=1; k=2; break;}
            case 1: { i=0; j=2; k=1; break;}
            case 2: { i=1; j=2; k=0; break;}
        }
        double alpha,beta;

        alpha = ( -D12[i] * D31[j] + D31[i] * D12[j] ) / D[index];
        beta = ( D42[i] * D12[j] - D12[i] * D42[j] ) / D[index];

        intersectingPoint[0]= x1[0] + (x3[0] -x1[0]) * beta;
        intersectingPoint[1]= x1[1] + (x3[1] -x1[1]) * beta;
        intersectingPoint[2]= x1[2] + (x3[2] -x1[2]) * beta;

        double EPSILON;
        EPSILON = std::numeric_limits<double>::epsilon();
        if (abs(D42[k] * alpha - D31[k] * beta - D12[k]) > 50 * EPSILON){
            std::cout << "FaceData - ERROR! The coordinates of the element do not lie on a plane, "<<abs(D42[k] * alpha - D31[k] * beta - D12[k]) << " is > "<<  50 * EPSILON <<" \n";
        }
        return intersectingPoint;
    }

    il::int_t FaceData::maxAbsArray(il::StaticArray<double, 3> x){
        // used to compute the centroid of a general quadrilateral
        // inputs
        //   -x array with 3 components
        // output
        //   -location: 0,1,2

        il::int_t index;
        for (il::int_t i = 0; i < 3; i++) {x[i] = abs(x[i]);}
        if (x[0] > x[1] ) {index = 0;} else {index = 1;}
        if (x[index] > x[2] ) {} else {index = 2;}
        return index;
    }
}