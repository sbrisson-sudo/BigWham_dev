//
// This file is part of BigWham.
//
// Created by Sylvain Brisson on 07/11/2025
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_POLYHEDRAL_TET_H
#define BIGWHAM_POLYHEDRAL_TET_H

#include "polyhedral.hh"

namespace bigwham {

template <int p> class Tetrahedron : public Polyhedral<p> {
public:
    Tetrahedron() : Polyhedral<p>() {

        // Vertices
        this->num_vertices_ = 4;
        this->vertices_.Resize(this->num_vertices_, 3);

        this->numFaces_ = 4;
        this->numVerticesPerFace_ = 3;

        // Assuming p=0
        this->num_nodes_ = 1;
        this->num_collocation_points_ = this->num_nodes_;

        this->collocation_points_.Resize(this->num_collocation_points_, 3);
        this->nodes_.Resize(this->num_nodes_, 3);

        // this->collocation_points_.Resize(this->num_collocation_points_, 3);
        // this->nodes_.Resize(this->num_nodes_, 3);
  }
  ~Tetrahedron() {}

  void setFaceIndices() override;
};

template <int p> 
void Tetrahedron<p>::setFaceIndices(){

    // Given as column-major
    il::Array2D<size_t> face_indices{il::value, {
        {0, 0, 0, 1},
        {1, 3, 2, 2},
        {3, 2, 1, 3},
    }};

    this->faceIndices_ = face_indices;
}


} // namespace bigwham


#endif // BIGWHAM_POLYHEDRAL_TET_H