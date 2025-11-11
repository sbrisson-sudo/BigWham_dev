//
// This file is part of BigWham.
//
// Created by Sylvain Brisson on 07/11/2025
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_POLYHEDRAL_HEX_H
#define BIGWHAM_POLYHEDRAL_HEX_H

#include "polyhedral.hh"

namespace bigwham {

template <int p> class Hexahedron : public Polyhedral<p> {
public:
    Hexahedron() : Polyhedral<p>() {

        // Vertices
        this->num_vertices_ = 8;
        this->vertices_.Resize(this->num_vertices_, 3);

        this->numFaces_ = 6;
        this->numVerticesPerFace_ = 4;

        // Assuming p=0
        this->num_nodes_ = 1;
        this->num_collocation_points_ = this->num_nodes_;

        this->collocation_points_.Resize(this->num_collocation_points_, 3);
        this->nodes_.Resize(this->num_nodes_, 3);

        // this->collocation_points_.Resize(this->num_collocation_points_, 3);
        // this->nodes_.Resize(this->num_nodes_, 3);
  }
  ~Hexahedron() {}

  void setFaceIndices() override;
};

template <int p> 
void Hexahedron<p>::setFaceIndices(){

    // Given as column-major
    il::Array2D<size_t> face_indices{il::value, {
        {1, 5, 4, 0, 0, 3},
        {2, 6, 7, 3, 1, 7},
        {6, 7, 3, 2, 5, 6},
        {5, 4, 0, 1, 4, 2}
    }};

    this->faceIndices_ = face_indices;
}


} // namespace bigwham


#endif // BIGWHAM_POLYHEDRAL_HEX_H