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

    // base class for boundary element
    template<int dim, int nvert, int p>
    class BoundaryElement {
    protected:

        int spatial_dimension_ = dim;           // spatial dimension
        int interpolation_order_ = p;           // order of interpolation for field on the element
        il::StaticArray2D<double, nvert, dim> vertices_;     // vertices' coordinates in global reference system -
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

        int getNumberOfVertices() const { return this->vertices_.size(0); };

        il::StaticArray<double, dim> getCentroid() const { return centroid_; };
        il::StaticArray<double, dim> getNormal() const { return n_; };
        il::StaticArray<double, dim> getTangent_1() const { return s_; };
        il::StaticArray<double, dim> getTangent_2() const { return t_; };
        il::StaticArray2D<double, dim, dim> rotationMatrix() const {
            il::StaticArray2D<double, dim, dim> R_{0.};
        };

        il::Array2D<double> getCollocationPoints() const { return collocation_points_; };
        virtual il::int_t getNumberOfCollocationPoints() const { return collocation_points_.size(0); };
        il::Array2D<double> getNodes() const { return nodes_; };
        virtual il::int_t getNumberOfNodes() const { return nodes_.size(0); };
        il::int_t getSpatialDimension() const { return spatial_dimension_; };
        il::int_t getInterpolationOrder() const { return interpolation_order_; };
    };

    template<int dim, int nvert, int p>
    BoundaryElement<dim, nvert, p>::BoundaryElement() = default;

    template<int dim, int nvert, int p>
    BoundaryElement<dim, nvert, p>::~BoundaryElement() = default;

};


#endif //BIGWHAM_BOUNDARYELEMENT_H
