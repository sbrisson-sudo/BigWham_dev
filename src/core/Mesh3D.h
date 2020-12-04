//
// This file is part of BigWham.
//
// Created by Alexis Sáez on 08.06.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIE_MESH3D_H
#define BIE_MESH3D_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray2D.h>
#include <elasticity/3d/constants.h>

#include <src/core/TriangularElementData.h>

namespace bie {

    class Mesh3D { // class for 2D triangular elements

    private:

        il::Array2D<double> coordinates_; // vertices' coordinates - size: number of vertices x 3
        il::Array2D<il::int_t> connectivity_; // connectivity matrix - size: number of elements x  3
        il::int_t interpolation_order_; // interpolation order: 0, 1 or 2

    public:
        //////////////////////////////////////////////////////////////////////////
        //        CONSTRUCTOR
        //////////////////////////////////////////////////////////////////////////

        //   Mesh3D default constructor;
        Mesh3D(){};
        ~Mesh3D(){};

        // Basic constructor with coordinates and connectivity matrices and interpolation order 'p'
        Mesh3D(const il::Array2D<double> &coordinates, const il::Array2D<il::int_t> &connectivity,
             const il::int_t interpolationOrder);

        //////////////////////////////////////////////////////////////////////////
        //        get-set functions  - i.e. public interfaces
        //////////////////////////////////////////////////////////////////////////

        il::int_t numberElts() const; // total number of elements
        il::int_t numberVertices() const; // total number of vertices
        il::int_t numberCollPtsElt() const; // number of collocation points per element
        il::int_t numberCollPts() const; // total number of collocation points = total number of nodes
        il::int_t interpolationOrder() const; // p: interpolation order

        ////////////////////////////////////////////////////////////////////////////////////////////
        //   Methods
        ////////////////////////////////////////////////////////////////////////////////////////////

        il::StaticArray2D<double,3,3> getVerticesElt(il::int_t ne) const; // vertices' coordinates for element
        // number/ID 'ne'
        bie::TriangularElementData getElementData(il::int_t ne) const; // build triangular element object for
        // element number/ID 'ne'
        il::Array2D<double> getCollocationPoints(); // get all the collocation points of the mesh
        il::Array2D<double> getNodes(); // get all the nodes of the mesh
        int getConnectivity(il::int_t e, il::int_t ln); // keep it here because it is passed to mma
    };

}
#endif //BIE_MESH3D_H

