//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_BOUNDARYELEMENT_H
#define BIGWHAM_BOUNDARYELEMENT_H

#include "il/container/1d/Array.h"
#include <limits>

#include <il/linearAlgebra.h>
#include <il/math.h>

#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linearAlgebra/dense/norm.h>

namespace bie {

// base class for boundary element
class BoundaryElement {

protected:
  using Real1D = il::Array<double>;
  using Real2D = il::Array2D<double>;
  int spatial_dimension_;   // spatial dimension
  int interpolation_order_; // order of interpolation for field on the element
  Real1D centroid_; // centroid of the element in global system of coordinates
  Real1D n_; // unit vector normal to element in global system of coordinates
  Real1D s_; // unit vector tangent to element in global system of coordinates,
  // in the direction from vertex 0 to vertex 1
  Real1D t_; // unit vector tangent to element in global system of coordinates,
  // orthogonal to s_ and n_ (un-used for 2D element)
  Real2D collocation_points_; // collocation points' coordinates in
                              // global reference system
  Real2D nodes_;              // nodes' coordinates in global reference system -
                              // size: number of nodes x dim
  Real2D vertices_; // vertices' coordinates in global reference system -
  Real2D rotMat_;   // rotation matrix for transformaion

  int n_vertices_;
  int n_nodes_;

public:
  BoundaryElement();
  BoundaryElement(int dim, int p) {
    this->centroid_.Resize(dim);
    this->n_.Resize(dim);
    this->s_.Resize(dim);
    this->t_.Resize(dim);
    this->spatial_dimension_ = dim;
    this->interpolation_order_ = p;
    this->rotMat_.Resize(dim, dim);
  };
  ~BoundaryElement();

  Real1D getCentroid() const { return this->centroid_; };
  Real1D getNormal() const { return this->n_; };
  Real1D getTangent_1() const { return this->s_; };
  Real1D getTangent_2() const { return this->t_; };

  Real2D getVertices() const { return this->vertices_; };
  Real2D getCollocationPoints() const { return this->collocation_points_; };
  Real2D getNodes() const { return nodes_; };

  il::int_t getNumberOfNodes() const { return this->nodes_.size(0); };
  il::int_t getSpatialDimension() const { return spatial_dimension_; };
  il::int_t getInterpolationOrder() const { return interpolation_order_; };
  il::int_t getNumberOfVertices() const { return this->vertices_.size(0); };
  il::int_t getNumberOfCollocationPoints() const {
    return collocation_points_.size(0);
  };

  virtual void setrotationMatrix() = 0;
  virtual void setElement(const Real2D &) = 0;
};

}; // namespace bie

#endif // BIGWHAM_BOUNDARYELEMENT_H
