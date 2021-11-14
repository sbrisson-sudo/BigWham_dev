//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#ifndef BIGWHAM_HMAT_H
#define BIGWHAM_HMAT_H

#pragma once

#include <hmat/hmatrix/LowRank.h>
#include <hmat/hmatrix/toHPattern.h>
#include <hmat/arrayFunctor/MatrixGenerator.h>

namespace bie {

template <il::int_t p, typename T>
class Hmat {
// this is a new Hmat class wich does not conntains the matrix-generator
// construction from the pattern built from the block cluster tree
// openmp parallel construction
// openmp parallel mat_vect dot product (non-permutted way)
 private:

  bie::HPattern pattern_; // the Hmat pattern

  il::int_t size0_; //n_row
  il::int_t size1_; //n_cols

  std::vector<bie::LowRank<T>> low_rank_blocks_; // vector of low rank blocks
  std::vector<il::Array2D<T>> full_rank_blocks_; // vector of full rank blocks

  bool built_= false;
  bool built_LR_=false;
  bool built_FR_=false;

  public:

   Hmat(const bie::HPattern& pattern){
     pattern_=pattern;
   };

   // -----------------------------------------------------------------------------
   void buildFR(const il::MatrixGenerator<T>& matrix_gen){
     // constructing the full rank blocks

  std::cout << " Loop on full blocks construction  \n";
  std::cout << " N full blocks "<< pattern_.n_FRB << " \n";

#pragma omp parallel
  {
    std::vector<il::Array2D<T>> private_full_rank_blocks;
#pragma omp for nowait schedule(static)
    for (il::int_t i = 0; i < pattern_.n_FRB; i++) {
      il::int_t i0 = pattern_.FRB_pattern(1, i);
      il::int_t j0 = pattern_.FRB_pattern(2, i);
      il::int_t iend = pattern_.FRB_pattern(3, i);
      il::int_t jend = pattern_.FRB_pattern(4, i);

      const il::int_t ni = p * (iend - i0);
      const il::int_t nj = p * (jend - j0);

      il::Array2D<T> sub{ni, nj};
      matrix_gen.set(i0, j0, il::io, sub.Edit());
      private_full_rank_blocks.push_back(sub);
      //full_rank_blocks_.push_back(sub);
    }
#ifdef _OPENMP
#pragma omp for schedule(static) ordered
    for(int i=0; i<omp_get_num_threads(); i++) {
#pragma omp ordered
      full_rank_blocks_.insert(full_rank_blocks_.end(), private_full_rank_blocks.begin(), private_full_rank_blocks.end());
    }
#else
    full_rank_blocks_.insert(full_rank_blocks_.end(), private_full_rank_blocks.begin(), private_full_rank_blocks.end());
#endif
  }
  built_FR_=true;

  }
////////////////////////////////////////////////////////////////////////////////
///
/// \param matrix_gen
/// \param epsilon
void buildLR(const il::MatrixGenerator<T>& matrix_gen,const double epsilon){
    // constructing the low rank blocks

    std::cout << " Loop on low rank blocks construction  \n";
    std::cout << " N low rank blocks "<< pattern_.n_LRB << " \n";

#pragma omp parallel
    {
      std::vector<bie::LowRank<T>> private_low_rank_blocks;
#pragma omp for nowait schedule(static)
      for (il::int_t i = 0; i < pattern_.n_LRB; i++) {
        il::int_t i0 = pattern_.LRB_pattern(1, i);
        il::int_t j0 = pattern_.LRB_pattern(2, i);
        il::int_t iend = pattern_.LRB_pattern(3, i);
        il::int_t jend = pattern_.LRB_pattern(4, i);

        il::Range range0{i0, iend}, range1{j0, jend};

        bie::LowRank<T> lra = bie::adaptiveCrossApproximation<p>(
            matrix_gen, range0, range1,
            epsilon);  // we need a LRA generator similar to the Matrix generator...

        private_low_rank_blocks.push_back(lra);  // not thread safe
        // store the rank in the low_rank pattern
        pattern_.LRB_pattern(5, i) = lra.A.size(1);
      }
#ifdef _OPENMP
#pragma omp for schedule(static) ordered
      for(int i=0; i<omp_get_num_threads(); i++) {
#pragma omp ordered
        low_rank_blocks_.insert(low_rank_blocks_.end(), private_low_rank_blocks.begin(), private_low_rank_blocks.end());
      }
#else
      low_rank_blocks_.insert(low_rank_blocks_.end(), private_low_rank_blocks.begin(), private_low_rank_blocks.end());
#endif
    }
    built_LR_=true;
  }
  //-----------------------------------------------------------------------------
  // filling up the h-matrix sub-blocks
  void build(const il::MatrixGenerator<T>& matrix_gen,const double epsilon){
    IL_EXPECT_FAST(matrix_gen.blockSize()==p);
    size0_=matrix_gen.size(0);
    size1_=matrix_gen.size(1);

    buildFR(matrix_gen);
    buildLR(matrix_gen,epsilon);
    built_=built_FR_ && built_LR_;
  }
//-----------------------------------------------------------------------------
  // getting the nb of entries of the hmatrix
  il::int_t nbOfEntries(){
    IL_EXPECT_FAST(built_);
    il::int_t n=0;

    for (il::int_t i=0;i<pattern_.n_FRB;i++){
      il::Array2DView<double> a = full_rank_blocks_[i].view();
      n+=a.size(0)*a.size(1);
    }
    for (il::int_t i=0;i<pattern_.n_LRB;i++) {
      il::Array2DView<double> a = low_rank_blocks_[i].A.view();
      il::Array2DView<double> b = low_rank_blocks_[i].B.view();
      n+=a.size(0)*a.size(1)+b.size(0)*b.size(1);
    }
    return n;
  }
  //-----------------------------------------------------------------------------
 // getting the compression ratio of the hmatrix (double which is <=1)
  double compressionRatio(){
    auto nb_elts = static_cast<double>(nbOfEntries());
    return nb_elts / static_cast<double>(size0_*size1_);
  }

  //--------------------------------------------------------------------------
  // H-Matrix vector multiplication without permutation
  il::Array<T> matvec(const il::Array<T> & x){
    IL_EXPECT_FAST(x.size()==size1_);

    il::Array<T> y{size0_,0.};

  // loop on full rank
#pragma omp parallel
    {
    il::Array<T> yprivate{size0_,0.};
#pragma omp for nowait schedule(static)
    for (il::int_t i = 0; i < pattern_.n_FRB; i++) {
      il::int_t i0=pattern_.FRB_pattern(1,i);
      il::int_t j0=pattern_.FRB_pattern(2,i);
      il::int_t iend=pattern_.FRB_pattern(3,i);
      il::int_t jend=pattern_.FRB_pattern(4,i);

      il::Array2DView<T> a = full_rank_blocks_[i].view();
      il::ArrayView<T> xs = x.view(il::Range{j0*p, jend*p});
      il::ArrayEdit<T> ys = yprivate.Edit(il::Range{i0*p, iend*p});

      il::blas(1.0, a, xs, 1.0, il::io, ys);
  }
  // the reduction below may be improved ?
#pragma omp critical
    for (il::int_t j=0;j<y.size();j++){
      y[j]+=yprivate[j];
    }
  };

  /// loop on low rank
#pragma omp parallel
  {
    il::Array<T> yprivate{size0_,0.};
#pragma omp for nowait schedule(static)
    for (il::int_t i = 0; i < pattern_.n_LRB; i++) {
      il::int_t i0 = pattern_.LRB_pattern(1, i);
      il::int_t j0 = pattern_.LRB_pattern(2, i);
      il::int_t iend = pattern_.LRB_pattern(3, i);
      il::int_t jend = pattern_.LRB_pattern(4, i);

      il::Array2DView<double> a = low_rank_blocks_[i].A.view();
      il::Array2DView<double> b = low_rank_blocks_[i].B.view();

      il::ArrayView<double> xs = x.view(il::Range{j0 * p, jend * p});
      il::ArrayEdit<double> ys = yprivate.Edit(il::Range{i0 * p, iend * p});
      const il::int_t r = a.size(1);
      il::Array<double> tmp{r, 0.0};

      il::blas(1.0, b, il::Dot::None, xs, 0.0, il::io,
               tmp.Edit());  // Note here we have stored b (not b^T)
      il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
    }
      // the reduction below may be improved ?
#pragma omp critical
     for (il::int_t j=0;j<y.size();j++){
        y[j]+=yprivate[j];
      }
  }
  return y;
  }

};


}
#endif

#pragma clang diagnostic pop