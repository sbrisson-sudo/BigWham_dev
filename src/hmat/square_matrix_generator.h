//
// This file is part of BigWham.
//
// Created by Brice Lecampion on  12.01.2023
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//
#pragma once

#ifndef BIGWHAM_SQMATGENERATOR_H
#define BIGWHAM_SQMATGENERATOR_H

#include <string>

#include "core/bie_kernel.h"
#include "core/elastic_properties.h"
#include "core/mesh.h"
#include "hmat/hierarchical_representation.h"
#include "hmat/arrayFunctor/matrix_generator.h"

namespace bigwham {

/*

 Class Square Matrix generator for BIE - note that the element type of the
 source and receiver elements should be the same!

 El: Element Type, Triangle<0>
 Bie_def: Kernel Type, BieElastostatic<Tri0, Tri0, bigwham::ElasticKernelType::H>
*/
template <typename T> class SquareMatrixGenerator : public MatrixGenerator<T> {

private:
  const std::shared_ptr<Mesh> mesh_;
  const std::shared_ptr<BieKernel<T>> bie_kernel_;

  il::int_t dof_dimension_; // unknowns per nodes
  il::int_t size_;          // total square matrix of size_*size_
  il::int_t num_points_;    // size_ / block_size_

public:
  SquareMatrixGenerator(const std::shared_ptr<Mesh> &mesh,
                        const std::shared_ptr<BieKernel<T>> &bie_kernel,
                        const std::shared_ptr<HRepresentation> &hr);
  il::int_t size(il::int_t d) const override;
  il::int_t blockSize() const override;
  il::int_t sizeAsBlocks(il::int_t d) const override;
  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<T> M) const override;
};

template <typename T>
SquareMatrixGenerator<T>::SquareMatrixGenerator(
    const std::shared_ptr<Mesh> &mesh,
    const std::shared_ptr<BieKernel<T>> &bie_kernel, const std::shared_ptr<HRepresentation> &hr)
    : mesh_(mesh), bie_kernel_(bie_kernel) {
  this->hr_ = hr;
  num_points_ = this->mesh_->num_collocation_points();
  dof_dimension_ = this->bie_kernel_->dof_dimension();
  size_ = num_points_ * dof_dimension_;
};

template <typename T>
il::int_t SquareMatrixGenerator<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);
  return size_;
};

template <typename T> il::int_t SquareMatrixGenerator<T>::blockSize() const {
  return dof_dimension_;
}

template <typename T>
il::int_t SquareMatrixGenerator<T>::sizeAsBlocks(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);
  return num_points_;
}

template <typename T>
void SquareMatrixGenerator<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                                   il::Array2DEdit<T> M) const {
  IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
  IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
  IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= num_points_);
  IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= num_points_);

  il::int_t jj = M.size(1) / blockSize();
#pragma omp parallel if (M.size(1) / blockSize() > 200)
  {
#pragma omp for
    for (il::int_t j1 = 0; j1 < M.size(1) / blockSize(); ++j1) {

      // El source_elt = El();
      // El receiver_elt = El();
      il::int_t k1 = b1 + j1;
      // j1 source node
      // from k1 - permute back to original mesh ordering using permutation
      // of
      // the clusters.
      il::int_t old_k1 = this->hr_->permutation_0_[k1];
      il::int_t e_k1 = this->mesh_->GetElementId(old_k1);
      il::int_t is_l = this->mesh_->GetElementCollocationId(old_k1);

      auto source_element = mesh_->GetElement(e_k1);

      for (il::int_t j0 = 0; j0 < M.size(0) / blockSize(); ++j0) {
        il::int_t k0 = b0 + j0;
        il::int_t old_k0 = this->hr_->permutation_1_[k0];
        il::int_t e_k0 = mesh_->GetElementId(old_k0); //  receiver element
        il::int_t ir_l = mesh_->GetElementCollocationId(old_k0);

        auto receiver_element = this->mesh_->GetElement(e_k0);
        std::vector<double> st = this->bie_kernel_->influence(
            *source_element, is_l, *receiver_element,
            ir_l); // column major
        // std::cout << "kernel size =" << st.size() << std::endl;
        // std::cout << "DOF dimension =" << dof_dimension_ << std::endl;
        IL_EXPECT_FAST(st.size() == dof_dimension_ * dof_dimension_);
        il::int_t k = 0;
        for (il::int_t j = 0; j < dof_dimension_; j++) {
          for (il::int_t i = 0; i < dof_dimension_; i++) {
            M(j0 * dof_dimension_ + i, j1 * dof_dimension_ + j) = st[k];
            k++;
          }
        }
      }
    }
  }
}

} // namespace bigwham

#endif // BIGWHAM_SQMATGENERATOR_H
