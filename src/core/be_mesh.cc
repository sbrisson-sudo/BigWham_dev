//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 01.02.23.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include "core/be_mesh.h"
#include "elements/rectangle.h"
#include "elements/segment.h"
#include "elements/triangle.h"
#include "elements/point.h"

namespace bigwham {
/* -------------------------------------------------------------------------- */

template <class ElemType>
il::Array2D<double> BEMesh<ElemType>::vertices(il::int_t element_id) const {
  il::Array2D<double> vertices{this->num_vertices_, this->spatial_dimension_};
  // loop over the vertices
  for (il::int_t j = 0; j < spatial_dimension_; j++) {
    for (il::int_t i = 0; i < num_vertices_; i++) {
      vertices(i, j) = coordinates_(connectivity_(element_id, i), j);
    }
  }
  return vertices;
}
/* -------------------------------------------------------------------------- */

template <class ElemType> void BEMesh<ElemType>::ConstructMesh() {
  il::int_t j = 0;
  // std::cout << num_colloc_pts_per_element_ << "\n";
  // std::cout << num_elements_ << "\n";
  // std::cout << element_list_.size() << "\n";
  for (il::int_t i = 0; i < element_list_.size(); i++) {
    // define element
    std::shared_ptr<BoundaryElement> elem = std::make_shared<ElemType>();
    auto vertices = this->vertices(i);
    elem->SetElement(vertices);
    // std::cout << num_colloc_pts_per_element_ << "\n";
    this->element_list_[i] = std::move(elem);

    // std::cout << num_colloc_pts_per_element_ << "\n";
    // define collocation points
    auto elem_colloc_pts = this->element_list_[i]->collocation_points();
    // std::cout << num_colloc_pts_per_element_ << "\n";
    for (il::int_t k = 0; k < spatial_dimension_; k++) {
      for (il::int_t j1 = 0; j1 < num_colloc_pts_per_element_; j1++) {
        collocation_points_(j + j1, k) = elem_colloc_pts(j1, k);
      }
    }
    j = j + num_colloc_pts_per_element_;
  }
}
/* -------------------------------------------------------------------------- */

template class BEMesh<Triangle<0>>;
template class BEMesh<Triangle<1>>;
template class BEMesh<Triangle<2>>;
template class BEMesh<Segment<0>>;
template class BEMesh<Segment<1>>;
template class BEMesh<Rectangle<0>>;
template class BEMesh<Point<2>>;
template class BEMesh<Point<3>>;

} // namespace bigwham
