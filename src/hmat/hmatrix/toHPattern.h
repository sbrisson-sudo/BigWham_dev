//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_TOHPATTERN_H
#define BIGWHAM_TOHPATTERN_H

#pragma once

#include <il/Tree.h>
#include <hmat/hmatrix/HMatrix.h>

namespace bie{

// construct pattern
//
// functions to get total number of sub-matrix blocks... (independant of dof size)
il::int_t nblocks_rec(const il::Tree<bie::SubHMatrix, 4>& tree, il::spot_t st) {
  const bie::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case bie::HMatrixType::FullRank: {
      return 1;
    } break;
    case bie::HMatrixType::LowRank: {
      return 1;
    } break;
    case bie::HMatrixType::Hierarchical: {
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
      return nblocks_rec(tree, st00)+nblocks_rec(tree, st01)+nblocks_rec( tree, st10)+nblocks_rec(tree, st11);
    } break;
    default:
      IL_UNREACHABLE;
  }
}

il::int_t nbBlocks(const il::Tree<bie::SubHMatrix, 4>& tree) {
  return nblocks_rec(tree,tree.root());
}

// functions to get total number of sub-matrix blocks... (independant of dof size)
il::int_t nfullblocks_rec(const il::Tree<bie::SubHMatrix, 4>& tree, il::spot_t st) {
  const bie::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case bie::HMatrixType::FullRank: {
      return 1;
    } break;
    case bie::HMatrixType::LowRank: {
      return 0;
    } break;
    case bie::HMatrixType::Hierarchical: {
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
      return nfullblocks_rec(tree, st00)+nfullblocks_rec(tree, st01)+nfullblocks_rec( tree, st10)+nfullblocks_rec(tree, st11);
    } break;
    default:
      IL_UNREACHABLE;
  }
}

il::int_t nbFullBlocks(const il::Tree<bie::SubHMatrix, 4>& tree) {
  return nfullblocks_rec(tree,tree.root());
}

// function to get the pattern
void pattern_rec(const il::Tree<bie::SubHMatrix, 4>& tree, il::spot_t st,
                 il::io_t,il::Array2D<il::int_t>& fr_pattern, il::int_t& nc_fr,il::Array2D<il::int_t>& lr_pattern, il::int_t& nc_lr) {
  const bie::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case bie::HMatrixType::FullRank: {
      fr_pattern(0,nc_fr)=st.index;
      fr_pattern(1,nc_fr)=info.range0.begin;
      fr_pattern(2,nc_fr)=info.range1.begin;
      fr_pattern(3,nc_fr)=info.range0.end;
      fr_pattern(4,nc_fr)=info.range1.end;
      fr_pattern(5,nc_fr)=0; //for full block - store the size
      nc_fr=nc_fr+1;
      return;
    }
    case bie::HMatrixType::LowRank: {
      lr_pattern(0,nc_lr)=st.index;
      lr_pattern(1,nc_lr)=info.range0.begin;
      lr_pattern(2,nc_lr)=info.range1.begin;
      lr_pattern(3,nc_lr)=info.range0.end;
      lr_pattern(4,nc_lr)=info.range1.end;
      lr_pattern(5,nc_lr)=1;
      nc_lr=nc_lr+1;
      return;
    }
    case bie::HMatrixType::Hierarchical: {
      const il::spot_t st00 = tree.child(st, 0);
      pattern_rec( tree, st00,  il::io, fr_pattern,nc_fr,lr_pattern,nc_lr);
      const il::spot_t st10 = tree.child(st, 1);
      pattern_rec( tree, st10,il::io, fr_pattern,nc_fr,lr_pattern,nc_lr);
      const il::spot_t st01 = tree.child(st, 2);
      pattern_rec( tree, st01, il::io, fr_pattern,nc_fr,lr_pattern,nc_lr);
      const il::spot_t st11 = tree.child(st, 3);
      pattern_rec( tree, st11, il::io, fr_pattern,nc_fr,lr_pattern,nc_lr);
      return;
    }
    default:
      IL_UNREACHABLE;
  }
}

// structure for storing the h-mat pattern - note irrespective of the number of dofs per nodes
struct HPattern {
  il::Array2D<il::int_t> FRB_pattern;   // full rank block pattern
  il::Array2D<il::int_t> LRB_pattern;  // low rank block pattern
  // pattern are stored as info a block k at pattern(0-5,k)
  //  spot,i_begin,j_begin,i_end,j_end,flag (0 for full rank, rank for low rank)
  il::int_t n_B{};   // number of blocks in the matrix pattern
  il::int_t n_FRB{}; // number of full rank blocks in the matrix pattern
  il::int_t n_LRB{}; // number of low rank blocks in the matrix pattern

  il::int_t nr{};  // total number of rows in the matrix pattern
  il::int_t nc{}; // total number of colums in the matrix pattern

};

// function to get the matrix pattern from the binary cluster tree
HPattern createPattern(const il::Tree<bie::SubHMatrix, 4>& tree){

  il::int_t nblocks = nbBlocks(tree);
  il::int_t n_full_blocks = nbFullBlocks(tree) ;
  il::int_t n_low_blocks = nblocks-n_full_blocks;

  HPattern my_pattern;
  my_pattern.n_B=nblocks;
  my_pattern.n_FRB =n_full_blocks;
  my_pattern.n_LRB=n_low_blocks;

  il::Array2D<il::int_t> fr_patt{6,n_full_blocks};
  il::Array2D<il::int_t> lr_patt{6,n_low_blocks};

  il::int_t nb_fr=0;
  il::int_t nb_lr=0;
  pattern_rec(tree,tree.root(),il::io,fr_patt,nb_fr,lr_patt,nb_lr);

  my_pattern.FRB_pattern=std::move(fr_patt);
  my_pattern.LRB_pattern=std::move(lr_patt);

  return my_pattern;
}


}

#endif  // BIGWHAM_TOHPATTERN_H
