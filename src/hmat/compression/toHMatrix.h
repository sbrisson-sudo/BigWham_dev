#pragma once

#include <il/Tree.h>

#include <hmat/arrayFunctor/MatrixGenerator.h>
#include <hmat/compression/adaptiveCrossApproximation.h>
#include <hmat/hmatrix/HMatrix.h>

#ifdef IL_PARALLEL
#include <tbb/tbb.h>
#endif

namespace il {

template <il::int_t p, typename T>
void hmatrix_rec(const il::MatrixGenerator<T>& matrix,
                 const il::Tree<il::SubHMatrix, 4>& tree, il::spot_t st,
                 double epsilon, il::spot_t shm, il::io_t, il::HMatrix<T>& hm) {
  const il::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case il::HMatrixType::FullRank: {
      const il::int_t n0 = p*(info.range0.end - info.range0.begin);
      const il::int_t n1 = p*(info.range1.end - info.range1.begin);
      hm.SetFullRank(shm, n0, n1);
      il::Array2DEdit<T> sub = hm.AsFullRank(shm);
      matrix.set(info.range0.begin, info.range1.begin, il::io, sub);
      return;
    } break;
    case il::HMatrixType::LowRank: {
      const il::int_t n0 = p*(info.range0.end - info.range0.begin);
      const il::int_t n1 = p*(info.range1.end - info.range1.begin);
      il::LowRank<T> sr = il::adaptiveCrossApproximation<p>(
          matrix, info.range0, info.range1, epsilon);
      const il::int_t r = sr.A.size(1);
      hm.SetLowRank(shm, n0, n1, r);
      il::Array2DEdit<T> A = hm.AsLowRankA(shm);
      for (il::int_t i1 = 0; i1 < A.size(1); ++i1) {
        for (il::int_t i0 = 0; i0 < A.size(0); ++i0) {
          A(i0, i1) = sr.A(i0, i1);
        }
      }
      il::Array2DEdit<T> B = hm.AsLowRankB(shm);
      for (il::int_t i1 = 0; i1 < B.size(1); ++i1) {
        for (il::int_t i0 = 0; i0 < B.size(0); ++i0) {
          B(i0, i1) = sr.B(i1, i0);
        }
      }
      return;
    } break;
    case il::HMatrixType::Hierarchical: {
      hm.SetHierarchical(shm);
//#ifdef IL_PARALLEL
// This cannot be done for the time being because allocating new nodes
// is not thread-safe
//      tbb::parallel_invoke(
//          [&] {
//            const il::spot_t st00 = tree.child(st, 0);
//            const il::spot_t shm00 = hm.child(shm, 0, 0);
//            hmatrix_rec<p>(matrix, tree, st00, epsilon, shm00, il::io, hm);
//          },
//
//          [&] {
//            const il::spot_t st10 = tree.child(st, 1);
//            const il::spot_t shm10 = hm.child(shm, 1, 0);
//            hmatrix_rec<p>(matrix, tree, st10, epsilon, shm10, il::io, hm);
//          },
//          [&] {
//            const il::spot_t st01 = tree.child(st, 2);
//            const il::spot_t shm01 = hm.child(shm, 0, 1);
//            hmatrix_rec<p>(matrix, tree, st01, epsilon, shm01, il::io, hm);
//          },
//          [&] {
//            const il::spot_t st11 = tree.child(st, 3);
//            const il::spot_t shm11 = hm.child(shm, 1, 1);
//            hmatrix_rec<p>(matrix, tree, st11, epsilon, shm11, il::io, hm);
//          });
//#else
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t shm00 = hm.child(shm, 0, 0);
      hmatrix_rec<p>(matrix, tree, st00, epsilon, shm00, il::io, hm);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t shm10 = hm.child(shm, 1, 0);
      hmatrix_rec<p>(matrix, tree, st10, epsilon, shm10, il::io, hm);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t shm01 = hm.child(shm, 0, 1);
      hmatrix_rec<p>(matrix, tree, st01, epsilon, shm01, il::io, hm);
      const il::spot_t st11 = tree.child(st, 3);
      const il::spot_t shm11 = hm.child(shm, 1, 1);
      hmatrix_rec<p>(matrix, tree, st11, epsilon, shm11, il::io, hm);
//#endif
      return;
    } break;
    default:
      IL_UNREACHABLE;
  }
}  // namespace il

template <typename T>
il::HMatrix<T> toHMatrix(const il::MatrixGenerator<T>& matrix,
                         const il::Tree<il::SubHMatrix, 4>& tree,
                         double epsilon) {
  il::HMatrix<T> ans{};
  if (matrix.blockSize() == 1) {
    hmatrix_rec<1>(matrix, tree, tree.root(), epsilon, ans.root(), il::io, ans);
  } else if (matrix.blockSize() == 2) {
    hmatrix_rec<2>(matrix, tree, tree.root(), epsilon, ans.root(), il::io, ans);
  } else if (matrix.blockSize() == 3) {
      hmatrix_rec<3>(matrix, tree, tree.root(), epsilon, ans.root(), il::io, ans); // needed for 3D
  } else {
    IL_UNREACHABLE;
  }
  return ans;
}

// construct pattern
//
// functions to get total number of sub-matrix blocks... (independant of dof size)
il::int_t nblocks_rec(const il::Tree<il::SubHMatrix, 4>& tree, il::spot_t st) {
  const il::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case il::HMatrixType::FullRank: {
      return 1;
    } break;
    case il::HMatrixType::LowRank: {
      return 1;
    } break;
    case il::HMatrixType::Hierarchical: {
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
     return nblocks_rec(tree, st00)+nblocks_rec(tree, st01)+nblocks_rec( tree, st10)+nblocks_rec(tree, st11);
    } break;
    default:
      IL_UNREACHABLE;
  }
}  // namespace il

il::int_t nbBlocks(const il::Tree<il::SubHMatrix, 4>& tree) {
  return nblocks_rec(tree,tree.root());
}
// functions to get total number of sub-matrix blocks... (independant of dof size)
il::int_t nfullblocks_rec(const il::Tree<il::SubHMatrix, 4>& tree, il::spot_t st) {
  const il::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case il::HMatrixType::FullRank: {
      return 1;
    } break;
    case il::HMatrixType::LowRank: {
      return 0;
    } break;
    case il::HMatrixType::Hierarchical: {
      const il::spot_t st00 = tree.child(st, 0);
      const il::spot_t st10 = tree.child(st, 1);
      const il::spot_t st01 = tree.child(st, 2);
      const il::spot_t st11 = tree.child(st, 3);
      return nfullblocks_rec(tree, st00)+nfullblocks_rec(tree, st01)+nfullblocks_rec( tree, st10)+nfullblocks_rec(tree, st11);
    } break;
    default:
      IL_UNREACHABLE;
  }
}  // namespace il

il::int_t nbFullBlocks(const il::Tree<il::SubHMatrix, 4>& tree) {
  return nfullblocks_rec(tree,tree.root());
}


// function to get the pattern now to avoid the recursion and pave the way to a parallel construction
void pattern_rec(const il::Tree<il::SubHMatrix, 4>& tree, il::spot_t st,
                    il::io_t,il::Array2D<il::int_t>& fr_pattern, il::int_t& nc_fr,il::Array2D<il::int_t>& lr_pattern, il::int_t& nc_lr) {
  const il::SubHMatrix info = tree.value(st);
  switch (info.type) {
    case il::HMatrixType::FullRank: {
      fr_pattern(0,nc_fr)=st.index;
      fr_pattern(1,nc_fr)=info.range0.begin;
      fr_pattern(2,nc_fr)=info.range1.begin;
      fr_pattern(3,nc_fr)=info.range0.end;
      fr_pattern(4,nc_fr)=info.range1.end;
      fr_pattern(5,nc_fr)=0;
      nc_fr=nc_fr+1;
      return;
    } break;
    case il::HMatrixType::LowRank: {
      lr_pattern(0,nc_lr)=st.index;
      lr_pattern(1,nc_lr)=info.range0.begin;
      lr_pattern(2,nc_lr)=info.range1.begin;
      lr_pattern(3,nc_lr)=info.range0.end;
      lr_pattern(4,nc_lr)=info.range1.end;
      lr_pattern(5,nc_lr)=1;
      nc_lr=nc_lr+1;
      return;
    } break;
    case il::HMatrixType::Hierarchical: {
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
  //  spot,i_begin,j_begin,i_end,j_end,flag (0 for full rank, 1 for low rank)
  il::int_t n_B;   // number of blocks
  il::int_t n_FRB; // number of full rank blocks
  il::int_t n_LRB; // number of low rank blocks
};

// function to get the matrix pattern from the binary cluster tree
HPattern createPattern(const il::Tree<il::SubHMatrix, 4>& tree){

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



}  // namespace il