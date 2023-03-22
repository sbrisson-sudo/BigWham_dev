//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 13.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_BEMESH_H
#define BIGWHAM_BEMESH_H

#include <memory>
#include <vector>

#include "core/mesh.h"
#include "elements/boundary_element.h"

namespace bie {

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
    // std::cout << elem->get_number_collocation_points() << "\n";
    interpolation_order_ = elem->get_interpolation_order();
    num_vertices_ = elem->get_number_vertices();
    num_colloc_pts_per_element_ = elem->get_number_collocation_points();
  }
  // Basic constructor with  coordinates and connectivity array and
  // element type !
  BEMesh(const il::Array2D<double> &coordinates,
         const il::Array2D<il::int_t> &connectivity)
    : Mesh(coordinates, connectivity) {
    // *this = BEMesh();
    // std::cout << num_colloc_pts_per_element_ << "\n";
    std::shared_ptr<BoundaryElement> elem = std::make_shared<ElemType>();
    interpolation_order_ = elem->get_interpolation_order();
    num_vertices_ = elem->get_number_vertices();
    num_colloc_pts_per_element_ = elem->get_number_collocation_points();
    this->num_collocation_points_ =
        this->num_elements_ * num_colloc_pts_per_element_;
    // std::cout << num_collocation_points_ << "\n";
    // std::cout << num_elements_ << "\n";
    collocation_points_.Resize(num_collocation_points_, spatial_dimension_);
  }

  virtual il::int_t get_element_id(il::int_t matrix_index) const override {
    return il::floor(matrix_index / (num_colloc_pts_per_element_));
  }

  virtual il::int_t
  get_element_collocation_id(il::int_t matrix_index) const override {
    return (matrix_index % (num_colloc_pts_per_element_));
  }

  il::Array2D<double> get_vertices(il::int_t element_id) const;

  virtual void construct_mesh() override;
};

} // namespace bie

#endif // BIGWHAM_BEMESH_H

// il::int_t connectivity(il::int_t e, il::int_t i)
//     const { // element e, local coordinates i - return global nodes
//   return connectivity_(e, i);
// }

// for (il::int_t i = 0; i < num_elements_; i++) {
//   BoundaryElement elem = std::make_unique<ElemType>();
//     }

// IL_EXPECT_FAST(spatial_dimension_ == 2 || spatial_dimension_ == 3);
// IL_EXPECT_FAST(connectivity.size(1) ==
// element_def_.getNumberOfVertices());

//  we do not check consistency of Connectivity here - such that actually
//  there can be no connectivity (can be used as dummy mesh)
//    for list of observation points for example
// this->number_vertex_ = element_def_.getNumberOfVertices();
// interpolation_order_ = element_def_.getInterpolationOrder();
// nodes_per_element_ = element_def_.getNumberOfNodes();
// }

//////////////////////////////////////////////////////////////////////////
//        get-set functions  - i.e. public interfaces
//////////////////////////////////////////////////////////////////////////

// il::int_t interpolationOrder() const { return interpolation_order_; };
// il::int_t numberOfNodes() const { return nodes_per_element_; };
// il::int_t numberCollocationPoints() const {
//   return (element_def_.getNumberOfCollocationPoints()) * n_elts_;
// }

// get the connectivity of an element -> A StaticArray of size 2 here !
// il::StaticArray<il::int_t, 2> connectivity(il::int_t k) // const
// {
//   il::StaticArray<il::int_t, 2> temp;
//   for (il::int_t i = 0; i < connectivity_.size(1); i++) {
//     temp[i] = connectivity_(k, i);
//   }
//   return temp;
// };
//

// il::int_t numberDDDofsPerElt() const {
//   return nodes_per_element_ * spatial_dimension_;
// }

// il::int_t numberDDDofs() const {
//   return (numberOfElts() * nodes_per_element_ * spatial_dimension_);
// }

////////////////////////////////////////////////////////////////////////////////////////////
//   Methods
////////////////////////////////////////////////////////////////////////////////////////////

// il::Array2D<double> getVertices(il::int_t ne) const {
//   il::Array2D<double> vertElt{number_vertex_, spatial_dimension_};
//   // loop over the vertices
//   for (il::int_t j = 0; j < spatial_dimension_; j++) {
//     for (il::int_t i = 0; i < number_vertex_; i++) {
//       vertElt(i, j) = coordinates_(connectivity_(ne, i), j);
//     }
//   }
//   return vertElt;
// }

// void setCurrentElement(il::int_t ne) {
//   il::Array2D<double> xv{
//       number_vertex_,
//       spatial_dimension_,
//       0,
//   };
//   xv = this->getVertices(ne);
//   this->element_def_.set_element(xv);
// }
