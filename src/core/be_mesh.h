//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BEMESH_H
#define BIGWHAM_BEMESH_H

#include <memory>
#include <vector>

#include "core/mesh.h"
#include "elements/boundary_element.h"

namespace bigwham {

// Class for Boundary element mesh
// a BE mesh has a element with dimension -1 compared to the spatial
// dimension... the coordinates of the element vertex have dimension spatial
// dimension this should replace the 2D and 3D class into a single class....
template <class ElemType> // E is the element type
class BEMesh : public Mesh {
private:
  // number of vertex per element
  il::int_t num_vertices_;
  // Interpolation order
  il::int_t interpolation_order_;
  //  collocation points oer element
  il::int_t num_colloc_pts_per_element_;


public:
  BEMesh() {
    std::shared_ptr<BoundaryElement> elem = std::make_shared<ElemType>();
    // std::cout << elem->num_collocation_points() << "\n";
    interpolation_order_ = elem->interpolation_order();
    num_vertices_ = elem->num_vertices();
    num_colloc_pts_per_element_ = elem->num_collocation_points();
  }
  // Basic constructor with  coordinates and connectivity array and
  // element type !
  BEMesh(const il::Array2D<double> &coordinates,
         const il::Array2D<il::int_t> &connectivity)
    : Mesh(coordinates, connectivity) {
    // *this = BEMesh();
    // std::cout << num_colloc_pts_per_element_ << "\n";
    std::shared_ptr<BoundaryElement> elem = std::make_shared<ElemType>();
    interpolation_order_ = elem->interpolation_order();
    num_vertices_ = elem->num_vertices();
    num_colloc_pts_per_element_ = elem->num_collocation_points();
    this->num_collocation_points_ =
        this->num_elements_ * num_colloc_pts_per_element_;
    // std::cout << num_collocation_points_ << "\n";
    // std::cout << num_elements_ << "\n";
    collocation_points_.Resize(num_collocation_points_, spatial_dimension_);
  }

  virtual il::int_t GetElementId(il::int_t matrix_index) const override {
    return il::floor(matrix_index / (num_colloc_pts_per_element_));
  }

  virtual il::int_t
  GetElementCollocationId(il::int_t matrix_index) const override {
    return (matrix_index % (num_colloc_pts_per_element_));
  }

  il::Array2D<double> vertices(il::int_t element_id) const;

  virtual void ConstructMesh() override;
};

} // namespace bigwham

#endif // BIGWHAM_BEMESH_H

