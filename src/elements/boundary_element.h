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

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/blas.h>

namespace bie {

// base class for boundary element
class BoundaryElement {

protected:
  int spatial_dimension_;   // spatial dimension
  int interpolation_order_; // order of interpolation for field on the element
  il::Array<double>
      centroid_; // centroid of the element in global system of coordinates
  il::Array<double>
      n_; // unit vector normal to element in global system of coordinates
  il::Array<double>
      s_; // unit vector tangent to element in global system of coordinates,
  // in the direction from vertex 0 to vertex 1
  il::Array<double>
      t_; // unit vector tangent to element in global system of coordinates,
  // orthogonal to s_ and n_ (un-used for 2D element)
  il::Array2D<double> collocation_points_; // collocation points' coordinates in
                                           // global reference system
  il::Array2D<double> nodes_; // nodes' coordinates in global reference system -
                              // size: number of nodes x dim
  il::Array2D<double>
      vertices_; // vertices' coordinates in global reference system -
  // rotation matrix for transformaion v_local = dot(rotation_matrix, v_global)
  il::Array2D<double> rotation_matrix_;
  // rotation matrix transformed v_global = dot(rotation_matrix_T, v_local)
  il::Array2D<double> rotation_matrix_T_;

  int num_vertices_;
  int num_nodes_;

public:
  BoundaryElement();
  BoundaryElement(int dim, int p) {
    this->centroid_.Resize(dim);
    this->n_.Resize(dim);
    this->s_.Resize(dim);
    this->t_.Resize(dim);
    this->spatial_dimension_ = dim;
    this->interpolation_order_ = p;
    this->rotation_matrix_.Resize(dim, dim);
    this->rotation_matrix_T_.Resize(dim, dim);
  }
  ~BoundaryElement() {}

  il::Array<double> get_centroid() const { return this->centroid_; };
  il::Array<double> get_normal() const { return this->n_; };
  il::Array<double> get_tangent_1() const { return this->s_; };
  il::Array<double> get_tangent_2() const { return this->t_; };

  il::Array2D<double> get_vertices() const { return this->vertices_; };
  il::Array2D<double> get_collocation_points() const {
    return this->collocation_points_;
  };
  il::Array2D<double> get_nodes() const { return this->nodes_; };

  il::int_t get_number_nodes() const { return this->num_nodes_; };
  il::int_t get_spatial_dimension() const { return this->spatial_dimension_; };
  il::int_t get_interpolation_order() const { return this->interpolation_order_; };
  il::int_t get_number_vertices() const { return this->vertices_.size(0); };
  il::int_t get_number_collocation_points() const {
    return this->collocation_points_.size(0);
  };

  il::Array<double> convert_to_local(const il::Array<double> global_vec) const {
    return il::dot(this->rotation_matrix_, global_vec);
  }
  il::Array<double> convert_to_global(const il::Array<double> local_vec) const {
    return il::dot(this->rotation_matrix_T_, local_vec);
  }

  // implement here rotation and rotation matrix transform here
  virtual void set_rotation_matrix() = 0;
  virtual void set_element(const il::Array2D<double> &coords_vertices) = 0;
};

}; // namespace bie

#endif // BIGWHAM_BOUNDARYELEMENT_H
