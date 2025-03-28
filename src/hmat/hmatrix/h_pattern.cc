//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include "h_pattern.h"

namespace bigwham {

// construct pattern
//
// functions to get total number of sub-matrix blocks... (independant of dof
// size)
il::int_t nblocks_rec(const il::Tree<bigwham::SubHMatrix, 4> &tree, il::spot_t st) {
  const bigwham::SubHMatrix info = tree.value(st);
  switch (info.type) {
  case bigwham::HMatrixType::FullRank: {
    return 1;
  } break;
  case bigwham::HMatrixType::LowRank: {
    return 1;
  } break;
  case bigwham::HMatrixType::Hierarchical: {
    if (tree.hasChild(st, 3)){ // Block has four children
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
      return nblocks_rec(tree, st00) + nblocks_rec(tree, st01) +
            nblocks_rec(tree, st10) + nblocks_rec(tree, st11);
    } else { // Block has only two children
      const il::spot_t st0 = tree.child(st, 0);
      const il::spot_t st1 = tree.child(st, 1);
      return nblocks_rec(tree, st0) + nblocks_rec(tree, st1);
    }
    
  } break;
  default:
    IL_UNREACHABLE;
  }
}

il::int_t nbBlocks(const il::Tree<bigwham::SubHMatrix, 4> &tree) {
  return nblocks_rec(tree, tree.root());
}

// functions to get total number of sub-matrix blocks... (independant of dof
// size)
il::int_t nfullblocks_rec(const il::Tree<bigwham::SubHMatrix, 4> &tree,
                          il::spot_t st) {
  const bigwham::SubHMatrix info = tree.value(st);
  switch (info.type) {
  case bigwham::HMatrixType::FullRank: {
    return 1;
  } break;
  case bigwham::HMatrixType::LowRank: {
    return 0;
  } break;
  case bigwham::HMatrixType::Hierarchical: {
    if (tree.hasChild(st, 3)){ // Block has four children
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
      return nfullblocks_rec(tree, st00) + nfullblocks_rec(tree, st01) +
            nfullblocks_rec(tree, st10) + nfullblocks_rec(tree, st11);
    } else { // Block has only two children
      const il::spot_t st0 = tree.child(st, 0);
      const il::spot_t st1 = tree.child(st, 1);
      return nfullblocks_rec(tree, st0) + nfullblocks_rec(tree, st1);
    }
  } break;
  default:
    IL_UNREACHABLE;
  }
}

il::int_t nbFullBlocks(const il::Tree<bigwham::SubHMatrix, 4> &tree) {
  return nfullblocks_rec(tree, tree.root());
}

// recursive function to get the pattern
void pattern_rec(const il::Tree<bigwham::SubHMatrix, 4> &tree, il::spot_t st,
                 il::io_t, il::Array2D<il::int_t> &fr_pattern, il::int_t &nc_fr,
                 il::Array2D<il::int_t> &lr_pattern, il::int_t &nc_lr) {
  const bigwham::SubHMatrix info = tree.value(st);
  switch (info.type) {
  case bigwham::HMatrixType::FullRank: {

    // std::cout << "FR block of size " << info.range0.end-info.range0.begin << " x " << info.range1.end-info.range1.begin << std::endl;

    fr_pattern(0, nc_fr) = st.index;
    fr_pattern(1, nc_fr) = info.range0.begin;
    fr_pattern(2, nc_fr) = info.range1.begin;
    fr_pattern(3, nc_fr) = info.range0.end;
    fr_pattern(4, nc_fr) = info.range1.end;
    fr_pattern(5, nc_fr) = 0; // for full block - store the size
    nc_fr = nc_fr + 1;
    return;
  }
  case bigwham::HMatrixType::LowRank: {

    // std::cout << "LR block of size " << info.range0.end-info.range0.begin << " x " << info.range1.end-info.range1.begin << std::endl;

    lr_pattern(0, nc_lr) = st.index;
    lr_pattern(1, nc_lr) = info.range0.begin;
    lr_pattern(2, nc_lr) = info.range1.begin;
    lr_pattern(3, nc_lr) = info.range0.end;
    lr_pattern(4, nc_lr) = info.range1.end;
    lr_pattern(5, nc_lr) = 1;
    nc_lr = nc_lr + 1;
    return;
  }
  case bigwham::HMatrixType::Hierarchical: {

    if (tree.hasChild(st, 3)){ // Block has four children
        const il::spot_t st00 = tree.child(st, 0);
        pattern_rec(tree, st00, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        const il::spot_t st10 = tree.child(st, 1);
        pattern_rec(tree, st10, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        const il::spot_t st01 = tree.child(st, 2);
        pattern_rec(tree, st01, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        const il::spot_t st11 = tree.child(st, 3);
        pattern_rec(tree, st11, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        return;
    } else { // Block has only two children
        const il::spot_t st00 = tree.child(st, 0);
        pattern_rec(tree, st00, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        const il::spot_t st10 = tree.child(st, 1);
        pattern_rec(tree, st10, il::io, fr_pattern, nc_fr, lr_pattern, nc_lr);
        return;
    }
  }
  default:
    IL_UNREACHABLE;
  }
}

// function to get the matrix pattern from the binary cluster tree
HPattern createPattern(const il::Tree<bigwham::SubHMatrix, 4> &tree) {

  il::int_t nblocks = nbBlocks(tree);
  il::int_t n_full_blocks = nbFullBlocks(tree);
  il::int_t n_low_blocks = nblocks - n_full_blocks;

  HPattern my_pattern;
  my_pattern.n_B = nblocks;
  my_pattern.n_FRB = n_full_blocks;
  my_pattern.n_LRB = n_low_blocks;

  il::Array2D<il::int_t> fr_patt{6, n_full_blocks};
  il::Array2D<il::int_t> lr_patt{6, n_low_blocks};

  il::int_t nb_fr = 0;
  il::int_t nb_lr = 0;
  pattern_rec(tree, tree.root(), il::io, fr_patt, nb_fr, lr_patt, nb_lr);

  my_pattern.FRB_pattern = std::move(fr_patt);
  my_pattern.LRB_pattern = std::move(lr_patt);

  return my_pattern;
}

} // namespace bigwham
