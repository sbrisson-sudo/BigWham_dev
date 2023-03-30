//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications 5.2.21: Moving to std::unique_ptr (C. Peruzzo)

#if defined(__clang__) && !defined(FMT_ICC_VERSION)
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#endif

#ifndef BIGWHAM_HMAT_H
#define BIGWHAM_HMAT_H

#ifdef IL_OPENMP
#include <omp.h>
#endif

#include <vector>

#include <hmat/arrayFunctor/matrix_generator.h>
#include <hmat/compression/adaptiveCrossApproximation.h>
#include <hmat/hmatrix/LowRank.h>
#include <hmat/hmatrix/h_pattern.h>

#include "hmat/hierarchical_representation.h"

namespace bie {

template <typename T> class Hmat {
  // this is a new Hmat class wich does not contains the matrix-generator
  // (neither the mesh etc.) construction from the pattern built from the block
  // cluster tree openmp parallel construction openmp parallel mat_vect dot
  // product (non-permutted way)
public:
  void hmatMemFree();
  Hmat() {}
  Hmat(const bie::MatrixGenerator<T> &matrix_gen, const double epsilon_aca);
  ~Hmat() {}
  void toHmat(const MatrixGenerator<T> &matrix_gen, const double epsilon_aca);
  std::vector<T> diagonal();
  std::vector<T> diagonalOriginal(const il::Array<il::int_t> &permutation);
  double compressionRatio();
  bool isBuilt() const { return isBuilt_; }
  il::int_t size(int k) const { return size_[k]; }
  bie::HPattern pattern() { return hr_->pattern_; }
  il::int_t dofDimension() const { return dof_dimension_; }
  il::int_t nbOfEntries();
  void fullBlocksOriginal(const il::Array<il::int_t> &permutation, il::io_t,
                          il::Array<T> &val_list, il::Array<int> &pos_list);
  // H-Matrix vector multiplication without permutation
  il::Array<T> matvec(const il::Array<T> &x);
  std::vector<T> matvec(const std::vector<T> &x);
  std::vector<T> matvecOriginal(const il::Array<il::int_t> &permutation,
                                const std::vector<T> &x);

private:
  void build(const bie::MatrixGenerator<T> &matrix_gen, const double epsilon);
  void buildFR(const bie::MatrixGenerator<T> &matrix_gen);
  void buildLR(const bie::MatrixGenerator<T> &matrix_gen, const double epsilon);

private:
  std::shared_ptr<HRepresentation> hr_;

  // shall we store the permutation(s) in that object ?

  il::int_t dof_dimension_;            //  dof per collocation points
  il::StaticArray<il::int_t, 2> size_; // size of tot mat (row, cols)

  std::vector<std::unique_ptr<bie::LowRank<T>>>
      low_rank_blocks_; // vector of low rank blocks
  std::vector<std::unique_ptr<il::Array2D<T>>>
      full_rank_blocks_; // vector of full rank blocks

  bool isBuilt_ = false;
  bool isBuilt_LR_ = false;
  bool isBuilt_FR_ = false;
};

} // namespace bie

#endif // BIGWHAM_HMAT_H

#if defined(__clang__) && !defined(FMT_ICC_VERSION)
#pragma clang diagnostic pop
#endif
