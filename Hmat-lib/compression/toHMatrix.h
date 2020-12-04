#pragma once

#include <il/Tree.h>

#include <Hmat-lib/arrayFunctor/MatrixGenerator.h>
#include <Hmat-lib/compression/adaptiveCrossApproximation.h>
#include <Hmat-lib/hmatrix/HMatrix.h>

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

}  // namespace il