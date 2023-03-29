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

#include <hmat/hmatrix/LowRank.h>
#include <hmat/hmatrix/toHPattern.h>
#include <hmat/arrayFunctor/matrix_generator.h>
#include <hmat/compression/adaptiveCrossApproximation.h>

#include "hmat/hierarchical_representation.h"

namespace bie {

template <typename T> class Hmat {
  // this is a new Hmat class wich does not contains the matrix-generator
  // (neither the mesh etc.) construction from the pattern built from the block
  // cluster tree openmp parallel construction openmp parallel mat_vect dot
  // product (non-permutted way)
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

public:
  Hmat() = default;
  ~Hmat() = default;

  // delete the memory pointed by low_rank_blocks_ and  full_rank_blocks_
  void hmatMemFree() {
    for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
      this->full_rank_blocks_[i].reset();
    }
    for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
      this->low_rank_blocks_[i].reset();
    }
  };

  //   // simple constructor from pattern  -> todo remove
  //   Hmat(const bie::HPattern& pattern){
  //     pattern_=pattern;
  //   };

  // direct constructor
  Hmat(const bie::MatrixGenerator<T> &matrix_gen,
       double epsilon_aca) {
    // construction directly
    this->hr_ = matrix_gen.hr();
    il::Timer tt;
    tt.Start();
    this->build(matrix_gen, epsilon_aca);
    tt.Stop();
    std::cout << "Creation of hmat done in " << tt.time() << "\n";
    std::cout << "Compression ratio - " << this->compressionRatio() << "\n";
    std::cout << "Hmat object - built "
              << "\n";
  }

  // Main constructor
  void toHmat(const MatrixGenerator<T> &matrix_gen,
              const double epsilon_aca) {
    this->hr_ = matrix_gen.hr();
    // construction directly
    il::Timer tt;
    tt.Start();
    this->build(matrix_gen, epsilon_aca);
    tt.Stop();
    std::cout << "Creation of hmat done in " << tt.time() << "\n";
    std::cout << "Compression ratio - " << this->compressionRatio() << "\n";
    std::cout << "Hmat object - built "
              << "\n";
  }

  // -----------------------------------------------------------------------------
  void buildFR(const bie::MatrixGenerator<T> &matrix_gen) {
    // construction of the full rank blocks
    std::cout << "Loop on full blocks construction  \n";
    std::cout << "N full blocks " << hr_->pattern_.n_FRB << " \n";

#pragma omp parallel
    {
      std::vector<std::unique_ptr<il::Array2D<T>>> private_full_rank_blocks;

#pragma omp for nowait schedule(static)
      for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
        il::int_t i0 = hr_->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr_->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr_->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr_->pattern_.FRB_pattern(4, i);

        const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
        const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

        std::unique_ptr<il::Array2D<T>> a =
            std::make_unique<il::Array2D<T>>(ni, nj);
        matrix_gen.set(i0, j0, il::io, (*a).Edit());
        private_full_rank_blocks.push_back(std::move(a));
      }
#ifdef IL_OPENMP
#pragma omp for schedule(static) ordered
      for (int i = 0; i < omp_get_num_threads(); i++) {
#pragma omp ordered
        full_rank_blocks_.insert(
            full_rank_blocks_.end(),
            std::make_move_iterator(private_full_rank_blocks.begin()),
            std::make_move_iterator(private_full_rank_blocks.end()));
      }
#else
      full_rank_blocks_.insert(
          full_rank_blocks_.end(),
          std::make_move_iterator(private_full_rank_blocks.begin()),
          std::make_move_iterator(private_full_rank_blocks.end()));
#endif
    }
    isBuilt_FR_ = true;
  }
  ////////////////////////////////////////////////////////////////////////////////
  ///
  /// \param matrix_gen
  /// \param epsilon
  void buildLR(const bie::MatrixGenerator<T> &matrix_gen,
               const double epsilon) {
    // constructing the low rank blocks
    dof_dimension_ = matrix_gen.blockSize();
    std::cout << "Loop on low rank blocks construction  \n";
    std::cout << "N low rank blocks " << hr_->pattern_.n_LRB << " \n";
    std::cout << "dof_dimension : " << dof_dimension_ << "\n";
#pragma omp parallel
    {
      std::vector<std::unique_ptr<bie::LowRank<T>>> private_low_rank_blocks;

#pragma omp for nowait schedule(static)
      for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
        il::int_t i0 = hr_->pattern_.LRB_pattern(1, i);
        il::int_t j0 = hr_->pattern_.LRB_pattern(2, i);
        il::int_t iend = hr_->pattern_.LRB_pattern(3, i);
        il::int_t jend = hr_->pattern_.LRB_pattern(4, i);
        il::Range range0{i0, iend}, range1{j0, jend};

        // we need a LRA generator virtual template similar to the Matrix
        // generator... here we have an if condition for the LRA call dependent
        // on dof_dimension_
        bie::LowRank<T> lra;
        if (matrix_gen.blockSize() == 1) {
          lra = bie::adaptiveCrossApproximation<1>(matrix_gen, range0, range1,
                                                   epsilon);
        } else if (matrix_gen.blockSize() == 2) {
          lra = bie::adaptiveCrossApproximation<2>(matrix_gen, range0, range1,
                                                   epsilon);
        } else if (matrix_gen.blockSize() == 3) {
          lra = bie::adaptiveCrossApproximation<3>(matrix_gen, range0, range1,
                                                   epsilon);
        } else {
          IL_UNREACHABLE;
        }
        std::unique_ptr<bie::LowRank<T>> lra_p(
            new bie::LowRank<T>(std::move(lra)));

        // store the rank in the low_rank pattern
        hr_->pattern_.LRB_pattern(5, i) = (*lra_p).A.size(1);

        private_low_rank_blocks.push_back(
            std::move(lra_p)); // lra_p does not exist after such call
      }

#ifdef IL_OPENMP
#pragma omp for schedule(static) ordered
      for (int i = 0; i < omp_get_num_threads(); i++) {
#pragma omp ordered
        low_rank_blocks_.insert(
            low_rank_blocks_.end(),
            std::make_move_iterator(private_low_rank_blocks.begin()),
            std::make_move_iterator(private_low_rank_blocks.end()));
      }
#else
      low_rank_blocks_.insert(
          low_rank_blocks_.end(),
          std::make_move_iterator(private_low_rank_blocks.begin()),
          std::make_move_iterator(private_low_rank_blocks.end()));
#endif
    }
    isBuilt_LR_ = true;
  }
  //-----------------------------------------------------------------------------
  // filling up the h-matrix sub-blocks
  void build(const bie::MatrixGenerator<T> &matrix_gen, const double epsilon) {
    dof_dimension_ = matrix_gen.blockSize();
    size_[0] = matrix_gen.size(0);
    size_[1] = matrix_gen.size(1);
    buildFR(matrix_gen);
    buildLR(matrix_gen, epsilon);
    isBuilt_ = isBuilt_FR_ && isBuilt_LR_;
  }
  //-----------------------------------------------------------------------------
  bool isBuilt() const { return isBuilt_; };
  //
  il::int_t size(int k) const { return size_[k]; };
  //
  bie::HPattern pattern() { return hr_->pattern_; }; // returning the Hmat pattern

  il::int_t dofDimension() const { return dof_dimension_; };

  //-----------------------------------------------------------------------------
  // getting the nb of entries of the hmatrix
  il::int_t nbOfEntries() {
    IL_EXPECT_FAST(isBuilt_);
    il::int_t n = 0;
    for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
      il::Array2DView<double> a = (*full_rank_blocks_[i]).view();
      n += a.size(0) * a.size(1);
    }
    for (il::int_t i = 0; i < hr_->pattern_.n_LRB; i++) {
      il::Array2DView<double> a = (*low_rank_blocks_[i]).A.view();
      il::Array2DView<double> b = (*low_rank_blocks_[i]).B.view();
      n += a.size(0) * a.size(1) + b.size(0) * b.size(1);
    }
    return n;
  }
  //-----------------------------------------------------------------------------
  // getting the compression ratio of the hmatrix (double which is <=1)
  double compressionRatio() {
    auto nb_elts = static_cast<double>(nbOfEntries());
    return nb_elts / static_cast<double>(size_[0] * size_[1]);
  }

  //--------------------------------------------------------------------------
  // H-Matrix vector multiplication without permutation
  // in - il:Array<T>
  // out - il:Array<T>
  il::Array<T> matvec(const il::Array<T> &x) {
    IL_EXPECT_FAST(x.size() == size_[1]);
    il::Array<T> y{size_[0], 0.};
    il::int_t n_B = hr_->pattern_.n_FRB + hr_->pattern_.n_LRB;

#pragma omp parallel
    {
      il::Array<T> yprivate{size_[0], 0.};
#pragma omp for nowait schedule(static)
      for (il::int_t i = 0; i < n_B; i++) {
        if (i < hr_->pattern_.n_FRB) // loop on full rank
        {
          il::int_t i0 = hr_->pattern_.FRB_pattern(1, i);
          il::int_t j0 = hr_->pattern_.FRB_pattern(2, i);
          il::int_t iend = hr_->pattern_.FRB_pattern(3, i);
          il::int_t jend = hr_->pattern_.FRB_pattern(4, i);

          il::Array2DView<T> a = (*full_rank_blocks_[i]).view();
          il::ArrayView<T> xs =
              x.view(il::Range{j0 * dof_dimension_, jend * dof_dimension_});
          il::ArrayEdit<T> ys = yprivate.Edit(
              il::Range{i0 * dof_dimension_, iend * dof_dimension_});

          il::blas(1.0, a, xs, 1.0, il::io, ys);
        } else /// loop on low rank
        {
          il::int_t ii = i - hr_->pattern_.n_FRB;
          il::int_t i0 = hr_->pattern_.LRB_pattern(1, ii);
          il::int_t j0 = hr_->pattern_.LRB_pattern(2, ii);
          il::int_t iend = hr_->pattern_.LRB_pattern(3, ii);
          il::int_t jend = hr_->pattern_.LRB_pattern(4, ii);

          il::Array2DView<double> a = (*low_rank_blocks_[ii]).A.view();
          il::Array2DView<double> b = (*low_rank_blocks_[ii]).B.view();

          il::ArrayView<double> xs =
              x.view(il::Range{j0 * dof_dimension_, jend * dof_dimension_});
          il::ArrayEdit<double> ys = yprivate.Edit(
              il::Range{i0 * dof_dimension_, iend * dof_dimension_});
          const il::int_t r = a.size(1);
          il::Array<double> tmp{r, 0.0};

          il::blas(1.0, b, il::Dot::None, xs, 0.0, il::io,
                   tmp.Edit()); // Note here we have stored b (not b^T)
          il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
        }
      }
      // the reduction below may be improved ?
#pragma omp critical
      for (il::int_t j = 0; j < y.size(); j++) {
        y[j] += yprivate[j];
      }
    };
    return y;
  }
  //--------------------------------------------------------------------------
  // H-Matrix vector multiplication with permutation for rectangular matrix
  // cases (only 1 permutation) in & out as std::vector todo : write another one
  // for the case of 2 permutations (rect. mat cases (for source != receivers)
  std::vector<T> matvecOriginal(const il::Array<il::int_t> &permutation,
                                const std::vector<T> &x) {

    il::Array<T> z{static_cast<il::int_t>(x.size())};
    // permutation of the dofs according to the re-ordering sue to clustering
    il::int_t ncolpoints = this->size(1) / dof_dimension_;
    for (il::int_t i = 0; i < ncolpoints; i++) {
      for (int j = 0; j < dof_dimension_; j++) {
        z[dof_dimension_ * i + j] = x[dof_dimension_ * permutation[i] + j];
      }
    }
    il::Array<T> y = this->matvec(z);
    std::vector<T> yout;
    yout.assign(y.size(), 0.);
    // permut back
    for (il::int_t i = 0; i < ncolpoints; i++) {
      for (int j = 0; j < dof_dimension_; j++) {
        yout[dof_dimension_ * permutation[i] + j] = y[dof_dimension_ * i + j];
      }
    }
    return yout;
  }
  ////////////////////////////////////////////////
  // matvect in and outs as std::vector
  std::vector<T> matvec(const std::vector<T> &x) {
    il::Array<T> xil{static_cast<il::int_t>(x.size())};
    // todo find a better way to convert il::Array to std::vect and vice versa !
    for (long i = 0; i < xil.size(); i++) {
      xil[i] = x[i];
    }
    il::Array<T> yil = this->matvec(xil);
    std::vector<T> y;
    y.reserve(static_cast<long>(yil.size()));
    for (long i = 0; i < yil.size(); i++) {
      y.push_back(yil[i]);
    }
    return y;
  }

  /////-----------------------------------------------------------------
  std::vector<T> diagonal() {
    // return diagonal in permutted state....
    IL_EXPECT_FAST(isBuilt_FR_);
    il::int_t diag_size = il::max(size_[0], size_[1]);
    std::vector<T> diag(static_cast<long>(diag_size), 0.);

    for (il::int_t k = 0; k < hr_->pattern_.n_FRB; k++) {
      il::Array2D<double> aux = *full_rank_blocks_[k];
      il::int_t i0 = hr_->pattern_.FRB_pattern(1, k);
      il::int_t j0 = hr_->pattern_.FRB_pattern(2, k);
      il::int_t iend = hr_->pattern_.FRB_pattern(3, k);
      il::int_t jend = hr_->pattern_.FRB_pattern(4, k);
      // check if it intersects the diagonal
      //
      bool in_lower = (i0 > j0) && (i0 > jend) && (iend > j0) && (iend > jend);
      bool in_upper = (i0 < j0) && (i0 < jend) && (iend < j0) && (iend < jend);
      if ((!in_lower) && (!in_upper)) // this fb intersect the diagonal....
      {
        for (il::int_t j = 0; j < aux.size(1); j++) {
          for (il::int_t i = 0; i < aux.size(0); i++) {
            if ((i + dof_dimension_ * i0) ==
                (j + dof_dimension_ * j0)) { // on diagonal !
              diag[(i + dof_dimension_ * i0)] = aux(i, j);
            }
          }
        }
      }
    }
    return diag;
  }
  /////
  std::vector<T> diagonalOriginal(const il::Array<il::int_t> &permutation) {
    // return diagonal in original state....
    il::int_t diag_size = il::max(size_[0], size_[1]);
    il::int_t ncolpoints = diag_size / dof_dimension_;
    std::vector<T> diag_raw = this->diagonal();
    std::vector<T> diag(diag_raw.size(), 0.);
    // permut back
    for (il::int_t i = 0; i < ncolpoints; i++) {
      for (int j = 0; j < dof_dimension_; j++) {
        diag[dof_dimension_ * permutation[i] + j] =
            diag_raw[dof_dimension_ * i + j];
      }
    }
    return diag;
  }

  //--------------------------------------------------------------------------
  void fullBlocksOriginal(const il::Array<il::int_t> &permutation, il::io_t,
                          il::Array<T> &val_list, il::Array<int> &pos_list) {
    // return the full blocks in the permutted Original dof state
    // in the val_list and pos_list 1D arrays.
    IL_EXPECT_FAST(isBuilt_FR_);
    IL_EXPECT_FAST(permutation.size() * dof_dimension_ == size_[1]);
    //  compute the number of  entries in the whole full rank blocks
    int nbfentry = 0;
    for (il::int_t i = 0; i < hr_->pattern_.n_FRB; i++) {
      il::Array2DView<double> aux = ((*full_rank_blocks_[i]).view());
      nbfentry = nbfentry + static_cast<int>(aux.size(0) * aux.size(1));
    }

    // prepare outputs
    pos_list.Resize(nbfentry * 2);
    val_list.Resize(nbfentry);

    il::Array<int> permutDOF{dof_dimension_ * permutation.size(), 0};
    IL_EXPECT_FAST(permutDOF.size() == size_[0]);
    for (il::int_t i = 0; i < permutation.size(); i++) {
      for (il::int_t j = 0; j < dof_dimension_; j++) {
        permutDOF[i * dof_dimension_ + j] = permutation[i] * dof_dimension_ + j;
      }
    }
    // loop on full rank and get i,j and val
    int nr = 0;
    int npos = 0;
    for (il::int_t k = 0; k < hr_->pattern_.n_FRB; k++) {
      il::Array2D<double> aux = *full_rank_blocks_[k];
      il::int_t i0 = hr_->pattern_.FRB_pattern(1, k);
      il::int_t j0 = hr_->pattern_.FRB_pattern(2, k);
      il::int_t index = 0;
      for (il::int_t j = 0; j < aux.size(1); j++) {
        for (il::int_t i = 0; i < aux.size(0); i++) {
          pos_list[npos + 2 * index] = permutDOF[(i + dof_dimension_ * i0)];
          pos_list[npos + 2 * index + 1] = permutDOF[(j + dof_dimension_ * j0)];
          val_list[nr + index] = aux(i, j);
          index++;
        }
      }
      nr = nr + static_cast<int>(aux.size(0) * aux.size(1));
      npos = npos + static_cast<int>(2 * aux.size(0) * aux.size(1));
    }
    std::cout << "done Full Block: nval " << val_list.size() << " / "
              << pos_list.size() / 2 << " n^2 "
              << (this->size_[0]) * (this->size_[1]) << "\n";
  }
};

} // namespace bie
#endif

#if defined(__clang__) && !defined(FMT_ICC_VERSION)
#pragma clang diagnostic pop
#endif
