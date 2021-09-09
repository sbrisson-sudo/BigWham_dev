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

template <typename T>
class Hmat {
 private:

  bie::HPattern pattern_;


  std::vector<il::LowRank<T>*> low_rank_blocks_;
  std::vector<il::int_t > lr_ranks_;             // list of ranks of low-rank
  std::vector<il::Array2D<T> * > full_rank_blocks_;

  public:

   Hmat(const bie::HPattern& pattern){
     pattern_=pattern;
   };
   //template <typename T>

  void buildFR(const il::MatrixGenerator<T>& matrix_gen){
  // building full rank blocks

  std::cout << " loop on fb construction  \n";

  for (il::int_t i=0;i<pattern_.n_FRB;i++){
    il::int_t i0=pattern_.FRB_pattern(i,1)*matrix_gen.blockSize();
    il::int_t i1=pattern_.FRB_pattern(i,2)*matrix_gen.blockSize();
    il::int_t j0=pattern_.FRB_pattern(i,3)*matrix_gen.blockSize();
    il::int_t j1=pattern_.FRB_pattern(i,4)*matrix_gen.blockSize();
    il::int_t ni=i1-i0;
    il::int_t nj=j1-j0;

    il::Array2D<T> sub{ni,nj};
    matrix_gen.set(i0, j0, il::io, sub.Edit());
   // std::cout << " fb #" << i <<"\n";

    full_rank_blocks_.push_back(&sub);

  }


  }

  void buildLR(){
    // building low rank
    il::int_t j=0;
  }

};


}
#endif  // BIGWHAM_HMAT_H
