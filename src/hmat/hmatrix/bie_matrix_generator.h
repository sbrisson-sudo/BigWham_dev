//
// This file is part of BigWham.
//
// Created by Brice Lecampion on  29.March.2023
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef BIGWHAM_BIE_MATRIX_GENERATOR_H
#define BIGWHAM_BIE_MATRIX_GENERATOR_H

#include "hmat/arrayFunctor/matrix_generator.h"
#include "hmat/hierarchical_representation.h"

#include "core/bie_kernel.h"
#include "core/elastic_properties.h"
#include "core/mesh.h"

namespace bie {

template <typename T> class BieMatrixGenerator : public MatrixGenerator<T> {

private:
  const il::Array<il::int_t> permutation0_; // row permutation
  const il::Array<il::int_t> permutation1_; // col permutation
  const std::shared_ptr<Mesh> mesh_src_;
  const std::shared_ptr<Mesh> mesh_rec_;
  const std::shared_ptr<BieKernel<T>> bie_kernel_src_;
  const std::shared_ptr<BieKernel<T>> bie_kernel_rec_;

  il::int_t block_size_; // dimension of the kernel (3 for 3D elasticity problems)
  il::int_t size0_;      // total rows of matrix
  il::int_t size1_;      // total cols of matrix
  il::int_t num_row_points;    // size0_ / block_size_
  il::int_t num_col_points;    // size1_ / block_size_

public:
BieMatrixGenerator(const std::shared_ptr<Mesh> &mesh_src,
const std::shared_ptr<Mesh> &mesh_rec,
                   const std::shared_ptr<BieKernel<T>> &bie_kernel_src,
                   const std::shared_ptr<BieKernel<T>> &bie_kernel_rec,
                   const HRepresentation &hr);
  virtual il::int_t size(il::int_t d) const override;
  virtual il::int_t blockSize() const override;
  virtual il::int_t sizeAsBlocks(il::int_t d) const override;
  virtual void set(il::int_t b0, il::int_t b1, il::io_t,il::Array2DEdit<T> M) const override;
};

template <typename T> inline BieMatrixGenerator<T>::BieMatrixGenerator(
    const std::shared_ptr<Mesh> &mesh,
    const std::shared_ptr<BieKernel<T>> &bie_kernel, const HRepresentation &hr)
    : mesh_(mesh), bie_kernel_(bie_kernel), permutation_(hr.permutation_0_) {
  num_points_ = this->mesh_->num_collocation_points();
  dof_dimension_ = this->mesh_->dof_dimension();
  size_ = num_points_ * dof_dimension_;
}

template <typename T> inline il::int_t BieMatrixGenerator<T>::size(il::int_t d) const {
  return size_;
}

template <typename T> inline il::int_t BieMatrixGenerator<T>::blockSize() const {
  return dof_dimension_;
}

template <typename T> inline il::int_t BieMatrixGenerator<T>::sizeAsBlocks(il::int_t d) const {
  return num_points_;
}

template <typename T> inline void BieMatrixGenerator<T>::set(il::int_t b0, il::int_t b1, il::io_t, il::Array2DEdit<T> M) const {

}

} // namespace bie
#endif  // BIGWHAM_BIE_MATRIX_GENERATOR_H