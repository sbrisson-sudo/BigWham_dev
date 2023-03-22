//
// This file is part of BigWham.
//
// Created by Alexis Sáez on 08.06.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include "Mesh3D.h"

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>
#include <iostream>

namespace bie {

    //////////////////////////////////////////////////////////////////////////
    //        CONSTRUCTOR
    //////////////////////////////////////////////////////////////////////////

    Mesh3D::Mesh3D(const il::Array2D<double> &coordinates, const il::Array2D<il::int_t> &connectivity,const il::int_t p, bool verbose){

        // inputs
        //   -coordinates: vertices' coordinates in the global reference system - size: number of vertices x 3
        //   -connectivity: connectivity matrix - size: number of elements x  3
        //   - p: interpolation order. It can be:
        //        0, for constant DD elements
        //        1, for linear DD elements
        //        2, for quadratic DD elements
        // output
        //   it builds the class object by initializing the inputs as member variables

        this->coordinates_ = coordinates;
        this->connectivity_ = connectivity;
        this->interpolation_order_ = p;
        if (verbose) {
            // to be printed in the mma kernel shell
            std::cout << " Mesh2D features" << "\n";
            std::cout << "  Number of vertices : " << this->numberVertices() << "\n";
            std::cout << "  Number of elements : " << this->numberOfElts() << "\n";
            std::cout << "  Number of collocation points : " << this->numberCollPts() << "\n";
         //   std::cout << "  Number of unknowns : " << 3 * this->numberCollPts() << "\n";
        }
    }

    //////////////////////////////////////////////////////////////////////////
    //        get-set functions  - i.e. public interfaces
    //////////////////////////////////////////////////////////////////////////

    il::int_t Mesh3D::numberOfElts() const { // total number of elements
        return this->connectivity_.size(0);
    }

    il::int_t Mesh3D::numberVertices() const { // total number of vertices
        return this->coordinates_.size(0);
    }

    il::int_t Mesh3D::numberCollPtsElt() const { // number of collocation points per element = number of
        // nodes per element
        il::int_t ncollpts;
        switch (this->interpolation_order_) {
            case 0: {
                ncollpts = 1;
                break;
            }
            case 1: {
                ncollpts = 3;
                break;
            }
            case 2: {
                ncollpts = 6;
            }
        }
        return ncollpts;
    }

    il::int_t Mesh3D::numberCollPts() const { // total number of collocation points = total number of nodes
        return this->numberOfElts() * this->numberCollPtsElt();
    }

    il::int_t Mesh3D::interpolationOrder() const {
        return this->interpolation_order_;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    //   Methods
    ////////////////////////////////////////////////////////////////////////////////////////////

    bie::FaceData Mesh3D::getElementData(const il::int_t ne) const {
        // inputs
        //   -ne: element number/ID
        // output
        //   -face element object of FaceData class containing all its member variables
        //   that are initialized when calling the constructor

        il::Array2D<double> xv=Mesh3D::getVerticesElt(ne);
        return FaceData(xv, this->interpolation_order_);
    }

    il::Array2D<double> Mesh3D::get_collocation_points() {
        // inputs
        //   -none; it makes use of coordinates_, connectivity_, and interpolation_order_, when calling
        //   using getElementData function
        // output
        //   -collocation points' coordinates of the whole mesh in the global reference system
        //   the output is delivered in the following 2D array form:
        //     x0 y0 z0
        //     x1 y1 z1
        //     .  .  .
        //     .  .  .
        //     .  .  .
        //     xN yN zN
        //    where N+1 is the total number of collocation points

        il::Array2D<double> collPts{this -> numberOfElts() * this -> numberCollPtsElt(), 3, 0};

        il::Array2D<double> collPtsElt;

        il::int_t j = 0;
        // loop over the elements
        for (il::int_t i = 0; i < this -> numberOfElts(); i++) {
            bie::FaceData elt = this->getElementData(i);
            collPtsElt = elt.get_collocation_points();
            // loop over collocation points per element
            for (il::int_t j1 = 0; j1 < this -> numberCollPtsElt(); j1++) {
                collPts(j, 0) = collPtsElt(j1, 0);
                collPts(j, 1) = collPtsElt(j1, 1);
                collPts(j, 2) = collPtsElt(j1, 2);
                j++;
            }
        }
        return collPts;
    }

    il::Array2D<double> Mesh3D::getNodes() {
        // inputs
        //   -none; it makes use of coordinates_, connectivity_, and interpolation_order_, when calling
        //   using getElementData function
        // output
        //   -nodes' coordinates of the whole mesh in the global reference system
        //   the output is delivered in the following 2D array form:
        //     x0 y0 z0
        //     x1 y1 z1
        //     .  .  .
        //     .  .  .
        //     .  .  .
        //     xN yN zN
        //    where N+1 is the total number of nodes

        il::Array2D<double> nodes{this -> numberOfElts() * this -> numberCollPtsElt(), 3, 0};

        il::Array2D<double> nodesElt;

        il::int_t j = 0;
        // loop over the elements
        for (il::int_t i = 0; i < this -> numberOfElts(); i++) {
            bie::FaceData elt = this->getElementData(i);
            nodesElt = elt.getNodes();
            // loop over collocation points per element
            for (il::int_t j1 = 0; j1 < this -> numberCollPtsElt(); j1++) {
                nodes(j, 0) = nodesElt(j1, 0);
                nodes(j, 1) = nodesElt(j1, 1);
                nodes(j, 2) = nodesElt(j1, 2);
                j++;
            }
        }
        return nodes;
    }

    int Mesh3D::getConnectivity(il::int_t e, il::int_t ln) {
        return connectivity_(e, ln);
    }

    il::Array2D<double> Mesh3D::getVerticesElt(il::int_t ne) const {
        // inputs
        //   -ne: element number/ID
        // output
        //   -vertices' coordinates of the element 'ne' in the global reference system
        //   the output is delivered in the following 2D array form:
        //     x0 y0 z0
        //     x1 y1 z1
        //     x2 y2 z2
        //     ...
        //
        //    where 0, 1, 2, ... are the vertices of the element
        // Note: The output expected by the methods that compute the quadratic kernel in 3D is
        // transposed w.r.t. this form, so you will find a transpose operation after calling this function
        // in the kernel files

        il::Array2D<double> vertElt{connectivity_.size(1),3};
        // loop over the vertices
        for (il::int_t i = 0; i < connectivity_.size(1); i++) {
            vertElt(i,0) = this -> coordinates_(this -> connectivity_(ne,i), 0); // x coordinate
            vertElt(i,1) = this -> coordinates_(this -> connectivity_(ne,i), 1); // y coordinate
            vertElt(i,2) = this -> coordinates_(this -> connectivity_(ne,i), 2); // z coordinate
        }
    return vertElt;
    }
}
