//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE
// file for more details.
//

#ifndef BIGWHAM_BOUNDARYELEMENT_H
#define BIGWHAM_BOUNDARYELEMENT_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/blas.h>

namespace bigwham {

// base class for boundary element
class BoundaryElement {

protected:
    // spatial dimension
  int spatial_dimension_;
    // order of interpolation for field on the element
  int interpolation_order_;
    // centroid of the element in global system of coordinates
  il::Array<double> centroid_;
    // unit vector normal to element in global system of coordinates
  il::Array<double>   normal_;
    // unit vector tangent to element in global system of coordinates,
    // in the direction from vertex 0 to vertex 1
  il::Array<double>   tangent1_;
    // unit vector tangent to element in global system of coordinates,
    // orthogonal to tangent1_ and normal_ (un-used for 2D element)
  il::Array<double>   tangent2_;
    // collocation points' coordinates in
    // global reference system
  il::Array2D<double> collocation_points_;
    // nodes' coordinates in global reference system -
    // size: number of nodes x dim
  il::Array2D<double> nodes_;
    // vertices' coordinates in global reference system -
  il::Array2D<double>   vertices_;
    // rotation matrix for transformaion v_local = dot(rotation_matrix, v_global)
  il::Array2D<double> rotation_matrix_;
    // rotation matrix transformed v_global = dot(rotation_matrix_T, v_local)
  il::Array2D<double> rotation_matrix_t_;

  int num_vertices_;
  int num_nodes_;
  int num_collocation_points_;

  // length for 1d, area for 2d
  double size_;

public:
  BoundaryElement() {}
  BoundaryElement(const int dim, const int p) {
    this->centroid_.Resize(dim, 0.);
    this->normal_.Resize(dim, 0.);
    this->tangent1_.Resize(dim, 0.);
    this->tangent2_.Resize(dim, 0.);
    this->spatial_dimension_ = dim;
    this->interpolation_order_ = p;
    this->rotation_matrix_.Resize(dim, dim, 0.);
    this->rotation_matrix_t_.Resize(dim, dim, 0.);
  }
  ~BoundaryElement() {}

  il::Array<double> centroid() const { return this->centroid_; };
  il::Array<double> normal() const { return this->normal_; };
  il::Array<double> tangent1() const { return this->tangent1_; };
  il::Array<double> tangent2() const { return this->tangent2_; };
  double size() const { return size_; }

  il::Array2D<double> vertices() const { return this->vertices_; };
  il::Array2D<double> collocation_points() const {
    return this->collocation_points_;
  };
  il::Array2D<double> nodes() const { return this->nodes_; };

  il::int_t num_nodes() const { return this->num_nodes_; };
  il::int_t spatial_dimension() const { return this->spatial_dimension_; };
  il::int_t interpolation_order() const {
    return this->interpolation_order_;
  };
  il::int_t num_vertices() const { return this->num_vertices_; };
  il::int_t num_collocation_points() const {
    return this->num_collocation_points_;
  };

  il::Array2D<double> rotation_matrix() const { return rotation_matrix_; }
  il::Array2D<double> rotation_matrix_t() const { return rotation_matrix_t_; }

  virtual il::Array<double> ConvertToLocal(const il::Array<double>& global_vec) const {
    return il::dot(this->rotation_matrix_, global_vec);
  }

  virtual il::Array<double> ConvertToGlobal(const il::Array<double>& local_vec) const {
    return il::dot(this->rotation_matrix_t_, local_vec);
  }

  // implement here rotation and rotation matrix transform here
  virtual void SetRotationMatrices() = 0;
  virtual void SetElement(const il::Array2D<double> &coords_vertices) = 0;
};

}; // namespace bigwham

#endif // BIGWHAM_BOUNDARYELEMENT_H
