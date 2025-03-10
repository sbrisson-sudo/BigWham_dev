//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications 5.2.21: Moving to std::unique_ptr (C. Peruzzo)

#if defined(__clang__) && !defined(FMT_ICC_VERSION)
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#endif

#ifndef BIGWHAM_HMAT_H
#define BIGWHAM_HMAT_H

#include <vector>

#include "hmat/arrayFunctor/matrix_generator.h"
#include "hmat/compression/adaptiveCrossApproximation.h"
#include "hmat/hmatrix/LowRank.h"
#include "hmat/hmatrix/h_pattern.h"
#include "hmat/hierarchical_representation.h"

namespace bigwham {

template <typename T> class Hmat {
  // this is a new Hmat class wich does not contains the matrix-generator
  // (neither the mesh etc.) construction from the pattern built from the block
  // cluster tree openmp parallel construction openmp parallel mat_vect dot
  // product (non-permutted way)
private:
    void build(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon);
    void buildFR(const bigwham::MatrixGenerator<T> & matrix_gen);
    template <il::int_t dim> void buildLR(const bigwham::MatrixGenerator<T> & matrix_gen,const double epsilon);

    std::shared_ptr<HRepresentation> hr_;

    il::int_t dof_dimension_;            //  dof per collocation points
    il::StaticArray<il::int_t, 2> size_; // size of tot mat (row, cols)

    std::vector<std::unique_ptr<bigwham::LowRank<T>>> low_rank_blocks_; // vector of low rank blocks
    std::vector<std::unique_ptr<il::Array2D<T>>> full_rank_blocks_; // vector of full rank blocks

    bool isBuilt_ = false;
    bool isBuilt_LR_ = false;
    bool isBuilt_FR_ = false;

    int n_openMP_threads_=8;

//#if defined(BIGWHAM_OPENMP)
    int frb_chunk_size_{1};
    int lrb_chunk_size_{1};
//#endif

    bool verbose_ = true;

public:
  void hmatMemFree();
  Hmat() = default;
  Hmat(const bigwham::MatrixGenerator<T> &matrix_gen, const double epsilon_aca,const int n_openMP_threads=8, const bool verbose=true);
  Hmat(const std::string &filename);

  ~Hmat() = default;
  void toHmat(const MatrixGenerator<T> &matrix_gen, const double epsilon_aca);
  std::vector<T> diagonal();
  std::vector<T> diagonalOriginal();
  double compressionRatio();
  [[nodiscard]] bool isBuilt() const { return isBuilt_; }
  [[nodiscard]] il::int_t size(int k) const { return size_[k]; }
  bigwham::HPattern pattern() { return hr_->pattern_; }
  [[nodiscard]] il::int_t dofDimension() const { return dof_dimension_; }
  il::int_t nbOfEntries();
  void fullBlocksOriginal(il::io_t, il::Array<T> & val_list,il::Array<int> & pos_list);
  void fullBlocksPerm(il::io_t, il::Array<T> & val_list,il::Array<int> & pos_list);
  // H-Matrix vector multiplication without permutation
  // il::Array<T> matvec(const il::Array<T> &x);
  std::vector<T> matvec(const std::vector<T> & x);
  il::Array<T> matvec(il::ArrayView<T> x);
  std::vector<T> matvecOriginal(const std::vector<T> & x);
  il::Array<T> matvecOriginal(il::ArrayView<T> x);

  void matvecOriginal(const il::ArrayView<T> x, il::Array<T>& yout);
  void matvecOriginal(const il::ArrayView<T> x, il::ArrayEdit<T> yout);

  void writeToFile(const std::string & filename);
  void readFromFile(const std::string & filename);

  std::shared_ptr<HRepresentation> getRepresentation() { return hr_; }

};

} // namespace bigwham

#endif // BIGWHAM_HMAT_H

#if defined(__clang__) && !defined(FMT_ICC_VERSION)
#pragma clang diagnostic pop
#endif
