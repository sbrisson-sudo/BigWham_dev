//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_HMAT_H
#define BIGWHAM_HMAT_H

#pragma once

#include <hmat/hmatrix/LowRank.h>
#include <hmat/hmatrix/toHPattern.h>
#include <hmat/arrayFunctor/MatrixGenerator.h>

namespace bie {

template <il::int_t p, typename T>
class Hmat {
 private:

  bie::HPattern pattern_;

  std::vector<il::LowRank<T>> low_rank_blocks_;
 // std::vector<il::int_t >           lr_ranks_;    // list of ranks of low-rank
  std::vector<il::Array2D<T> > full_rank_blocks_;

  bool built_= false;
  bool built_LR_=false;
  bool built_FR_=false;

  public:

   Hmat(const bie::HPattern& pattern){
     pattern_=pattern;
   };
   // template <typename T>

   // For simplicity for now, we do not OpenMP the loop on full blocks creation
   // (note the multithreading is done at the matrix_gen.set level)
  void buildFR(const il::MatrixGenerator<T>& matrix_gen){
  // building full rank blocks
  IL_EXPECT_FAST(matrix_gen.blockSize()==p);

  std::cout << " Loop on full blocks construction  \n";
  std::cout << " N full blocks "<< pattern_.n_FRB << " \n";

  for (il::int_t i=0;i<pattern_.n_FRB;i++){

    il::int_t i0=pattern_.FRB_pattern(1,i);
    il::int_t j0=pattern_.FRB_pattern(2,i);
    il::int_t iend=pattern_.FRB_pattern(3,i);
    il::int_t jend=pattern_.FRB_pattern(4,i);

    const il::int_t ni=p*(iend-i0);
    const il::int_t nj=p*(jend-j0);

    il::Array2D<T> sub{ni,nj};
    matrix_gen.set(i0, j0, il::io, sub.Edit());
    full_rank_blocks_.push_back(sub); // not thread safe ? - do a loop to reserve first ?
  }

  built_FR_=true;

  }
////////////////////////////////////////////////////////////////////////////////
///
/// \param matrix_gen
/// \param epsilon
void buildLR(const il::MatrixGenerator<T>& matrix_gen, const double epsilon){
    // building low rank
    IL_EXPECT_FAST(matrix_gen.blockSize()==p);

    std::cout << " Loop on low rank blocks construction  \n";
    std::cout << " N low rank blocks "<< pattern_.n_LRB << " \n";

    for (il::int_t i=0;i<pattern_.n_LRB;i++) {

      il::int_t i0=pattern_.LRB_pattern(1,i);
      il::int_t j0=pattern_.LRB_pattern(2,i);
      il::int_t iend=pattern_.LRB_pattern(3,i);
      il::int_t jend=pattern_.LRB_pattern(4,i);

      il::Range range0{i0,iend},range1{j0,jend};

      il::LowRank<T> lra = il::adaptiveCrossApproximation<p>(
          matrix_gen, range0, range1, epsilon);   // we need a LRA generator similar to the Matrix generator...

      low_rank_blocks_.push_back(lra); // not thread safe
      // store the rank in the low_rank pattern
      pattern_.LRB_pattern(5,i)=lra.A.size(1);
    }
    built_LR_=true;
  }

  void build(const il::MatrixGenerator<T>& matrix_gen, const double epsilon){
    buildFR(matrix_gen);
    buildLR(matrix_gen,epsilon);
    built_=built_FR_ & built_LR_;
  }
  //--------------------------------------------------------------------------
  //// MAtrix vector multiplication
  il::Array<T> matvec(const il::Array<T> & x){
    //  without permutation at this point
  il::int_t npts = x.size()/p;
  il::Array<T> y{x.size(),0.};

  // loop on full rank
#ifndef NUMBEROFTHREADS
#define NUMBEROFTHREADS 4
#endif
#pragma omp parallel for num_threads(NUMBEROFTHREADS)
  for (il::int_t i = 0; i < pattern_.n_FRB; i++) {

    il::int_t i0=pattern_.FRB_pattern(1,i);
    il::int_t j0=pattern_.FRB_pattern(2,i);
    il::int_t iend=pattern_.FRB_pattern(3,i);
    il::int_t jend=pattern_.FRB_pattern(4,i);

    il::Array2DView<double> a = full_rank_blocks_[i].view();

    il::ArrayView<double> xs = x.view(il::Range{j0*p, jend*p});
    il::ArrayEdit<double> ys = y.Edit(il::Range{i0*p, iend*p});

    il::blas(1.0, a, xs, 1.0, il::io, ys);
  }

  /// loop on low rank
#pragma omp parallel for num_threads(NUMBEROFTHREADS)
    for (il::int_t i = 0; i < pattern_.n_LRB; i++) {

    il::int_t i0=pattern_.LRB_pattern(1,i);
    il::int_t j0=pattern_.LRB_pattern(2,i);
    il::int_t iend=pattern_.LRB_pattern(3,i);
    il::int_t jend=pattern_.LRB_pattern(4,i);

    il::Array2DView<double> a =  low_rank_blocks_[i].A.view() ;
    il::Array2DView<double> b =  low_rank_blocks_[i].B.view() ;

    il::ArrayView<double> xs = x.view(il::Range{j0*p, jend*p});
    il::ArrayEdit<double> ys = y.Edit(il::Range{i0*p, iend*p});
    const il::int_t r = a.size(1);
    il::Array<double> tmp{r, 0.0};

    il::blas(1.0, b, il::Dot::None, xs, 0.0, il::io, tmp.Edit()); // Note here we have stored b (not b^T)
    il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
  }

  return y;

  }

};


}
#endif  // BIGWHAM_HMAT_H
