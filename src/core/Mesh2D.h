//
// HFP project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#ifndef BIE_MESH_H
#define BIE_MESH_H

// include std libraries
#include <algorithm>
#include <cmath>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra/dense/norm.h>

// Inclusion from bie
#include <src/core/SegmentData.h>

namespace bie {

///// 1D mesh class
class Mesh {  // class for   1D segment elements

 private:
  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> coordinates_;

  // Connectivity matrix - size: number of elements x  2
  il::Array2D<il::int_t> connectivity_;

  // Interpolation order
  il::int_t interpolation_order_;

  // Dof handle matrices
  // for displacements  Discontinuities - size: number of elements x 2 dofs per
  // coordinates x (order interpolation + 1)
//  il::Array2D<il::int_t> dof_handle_dd_; //todo : get rid of this dof_handle

 public:
  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////////

  // todo: naming of the different entities are not consistent

  //   Mesh()default;
  Mesh(){};
  ~Mesh(){};

  // Basic constructor with  coordinates and connectivity array and
  // interpolation order
  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity,
       const il::int_t interpolationOrder) {
    // check validity of inputs

    IL_EXPECT_FAST(Coordinates.size(0) > 1 && Coordinates.size(1) == 2);
    // P0 and P1 elements only for now
    IL_EXPECT_FAST(interpolationOrder == 0 || interpolationOrder == 1);
    // check connectivity and coordinates consistency ??? currently no ->
    // they should be properly compatible

    coordinates_ = Coordinates;
    connectivity_ = Connectivity;
    interpolation_order_ = interpolationOrder;

    il::int_t  nelts=connectivity_.size(0);
    il::int_t p=interpolationOrder;

    /// Discontinuous Polynomial DOF handles
//    il::Array2D<il::int_t> id_dd{nelts, 2 * (p + 1), 0};
//
//    for (il::int_t j = 0; j < 2 * (p + 1); j++) {
//      for (il::int_t i = 0; i < nelts; i++) {
//        id_dd(i, j) = i * 2 * (p + 1) + j;
//      }
//    }
//    dof_handle_dd_ = id_dd;  /// dof

  };

  //////////////////////////////////////////////////////////////////////////
  //        get-set functions  - i.e. public interfaces
  //////////////////////////////////////////////////////////////////////////

  il::int_t numberOfElts() const { return connectivity_.size(0); };

  il::int_t numberOfNodes() const { return coordinates_.size(0); };

  // nodal coordinates related.
  il::Array2D<double> coordinates() const { return coordinates_; };

  // Read a particular element of the coordinates coordinates
  double coordinates(il::int_t k, il::int_t i)  // const
  {
    return coordinates_(k, i);
  }

  il::StaticArray<double, 2> coordinates(il::int_t k)  // const
  {
    il::StaticArray<double, 2> temp;
    temp[0] = coordinates_(k, 0);
    temp[1] = coordinates_(k, 1);
    return temp;
  };

  // connectivity related
  il::Array2D<il::int_t> connectivity() const { return connectivity_; };
  // get the connectivity of an element -> A StaticArray of size 2 here !
  il::StaticArray<il::int_t, 2> connectivity(il::int_t k)  // const
  {
    il::StaticArray<il::int_t, 2> temp;
    for (il::int_t i = 0; i < connectivity_.size(1); i++) {
      temp[i] = connectivity_(k, i);
    }
    return temp;
  };

  //
  il::int_t connectivity(il::int_t e, il::int_t i)  // const
  {
    // element e, local coordinates i - return global nodes
    return connectivity_(e, i);
  }

  // interpolation order
  il::int_t interpolationOrder() const { return interpolation_order_; }

  // dofs related.....
  il::int_t numberDDDofsPerElt() const {
    return connectivity_.size(1) * (interpolation_order_ + 1);
  }

  il::int_t numberDDDofs() const {
    return (numberOfElts() * (interpolation_order_ + 1) * 2);
  }

  il::int_t numberCollocationPoints() const {
    return (numberOfElts() * (interpolation_order_ + 1) );
  }

//  il::int_t dofDD(il::int_t k, il::int_t i) // const
//  {
//    // coordinates k, dof i -> return global equation iD
//    return dof_handle_dd_(k, i);  // element , dof dim.
//  }

//  il::Array2D<double> getNodes()
//  {
//
//  }


  ////////////////////////////////////////////////////////////////////////////////////////////
  //   Methods
  ////////////////////////////////////////////////////////////////////////////////////////////

  il::Array2D<double> getCollocationPoints();

  bie::SegmentData getElementData(il::int_t ne)  const;
};

}
#endif  // BIE_MESH_H
