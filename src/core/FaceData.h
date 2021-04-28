//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.01.21 based on a work done by Alexis Sáez on 08.06.20.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIE_FACEDATA_H
#define BIE_FACEDATA_H

#pragma once

#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/Array2D.h>

namespace bie {

class FaceData {

private:

    il::int_t interpolation_order_;    // interpolation order: 0, 1 or 2
    il::StaticArray<double, 3> n_;     // unit vector normal to element in global system of coordinates
    il::StaticArray<double, 3> s_;     // unit vector tangent to element in global system of coordinates,
                                       // in the direction from vertex 0 to vertex 1
    il::StaticArray<double, 3> t_;     // unit vector tangent to element in global system of coordinates,
                                       // orthogonal to s_ and n_
    il::StaticArray<double, 3> xc_;    // centroid of the element in global system of coordinates
    il::Array2D<double> vertices_;     // vertices' coordinates in global reference system -
                                       // size: number of vertices x 3
    il::int_t NoV_;                    // Number of vertices
    il::Array2D<double> nodes_;        // nodes' coordinates in global reference system - size: number of nodes x 3
    il::Array2D<double> collocation_points_; // collocation points' coordinates in global reference system
                                             // - size: number of collocation points x 3

    const double beta_ = 1.5 * 0.166666666666667; // collocation points' location parameter for linear
    // elements (see documentation)
    // collocation points' location parameters for quadratic elements
    const double beta1_ = 0.35; // 1.5 * 0.091576213509771 related to nodes at vertices (see documentation)
    const double beta2_ = 0.35; // 1.5 * 0.10810301816807 related to middle-edge nodes (see documentation)
public:

    //////////////////////////////////////////////////////////////////////////
    //        CONSTRUCTOR
    //////////////////////////////////////////////////////////////////////////

    // Basic constructor with coordinates of vertices 'xv' and interpolation order 'p'
    FaceData(il::Array2D<double> xv, il::int_t p);

    //////////////////////////////////////////////////////////////////////////
    //        get-set functions  - i.e. public interfaces
    //////////////////////////////////////////////////////////////////////////

    il::Array2D<double> getNodes();
    double getNoV();
    il::Array<double> getNormal();
    il::Array<double> getS1();
    il::Array<double> getS2();
    il::Array2D<double> getCollocationPoints();
    il::Array2D<double> getVertices(); // this function is a bit silly because
    // the object is indeed constructed by the vertices as input, however is needed due to the way the
    // construction of the elasticity matrix is coded up for quadratic (p=2) elements. This function
    // should be deleted in the future
    double getBeta1() const; // this function is needed because in the construction of the
    // elasticity matrix the computation of the collocation points is duplicated (done by a previous
    // function of Dmitry). This function should be deleted in the future
    double getBeta2() const; // same as before

    // uncomment the following if needed ...
//    il::StaticArray<double, 3> getCentroid() { return xc_;};
//    il::StaticArray<double, 3> getNormal() { return n_;};
//    il::StaticArray<double, 3> getTangent1() { return s_;};
//    il::StaticArray<double, 3> getTangent2() { return t_;};
//    const double getBeta() { return beta_;};

    //////////////////////////////////////////////////////////////////////////
    //   Methods
    //////////////////////////////////////////////////////////////////////////

        il::StaticArray<double, 3> computeCentroid(il::Array2D<double> xv);
        il::Array2D<double> computeNodes();
        il::Array2D<double> computeCollocationPoints();
        il::Array2D<double> rotationMatrix(bool Transposed = false); // or quadratic elements the function written by Dmitry is used
private:
        il::int_t maxAbsArray(il::StaticArray<double, 3> x);
        il::StaticArray<double, 3> intersectLines(il::StaticArray<double, 3> x1, il::StaticArray<double, 3> x2,
                                                  il::StaticArray<double, 3> x3, il::StaticArray<double, 3> x4);
    };
}

#endif  // BIE_FACEDATA_H
