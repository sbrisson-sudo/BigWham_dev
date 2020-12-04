
#pragma once

#include <tbb/tbb.h>

#include <il/Array.h>
#include <il/linearAlgebra.h>


#include <Hmat-lib/hmatrix/HMatrix.h>

namespace il {

template <typename T>
void dot_rec(il::int_t parallelism, const il::HMatrix<T>& A, il::spot_t s, il::ArrayView<T> x,
             il::io_t, il::ArrayEdit<T> y) {
  IL_EXPECT_FAST(A.size(0, s) == y.size());
  IL_EXPECT_FAST(A.size(1, s) == x.size());

  if (A.isFullRank(s)) {
    il::Array2DView<T> a = A.asFullRank(s);
    il::blas(1.0, a, x, 1.0, il::io, y);
    return;
  } else if (A.isLowRank(s)) {
    il::Array2DView<T> a = A.asLowRankA(s);
    il::Array2DView<T> b = A.asLowRankB(s);
    const il::int_t r = a.size(1);
    il::Array<T> tmp{r, 0.0};
    il::blas(1.0, b, il::Dot::Transpose, x, 0.0, il::io, tmp.Edit());
    il::blas(1.0, a, tmp.view(), 1.0, il::io, y);
    return;
  } else if (A.isHierarchical(s)) {
    const il::spot_t s00 = A.child(s, 0, 0);
    const il::spot_t s10 = A.child(s, 1, 0);
    const il::spot_t s01 = A.child(s, 0, 1);
    const il::spot_t s11 = A.child(s, 1, 1);
    const il::int_t n00 = A.size(0, s00);
    const il::int_t n10 = A.size(0, s10);
    const il::int_t n01 = A.size(1, s00);
    const il::int_t n11 = A.size(1, s01);
    il::ArrayView<T> x0 = x.view(il::Range{0, n01});
    il::ArrayView<T> x1 = x.view(il::Range{n01, n01 + n11});
    il::ArrayEdit<T> y0 = y.Edit(il::Range{0, n00});
    il::ArrayEdit<T> y1 = y.Edit(il::Range{n00, n00 + n10});


    if (parallelism >= 0) {
      il::Array<T> y0_par{n00, 0.0};
      il::Array<T> y1_par{n10, 0.0};
      tbb::parallel_invoke(
          [&]() {
            dot_rec(parallelism - 1, A, s00, x0, il::io, y0);
          },
          [&]() {
            dot_rec(parallelism - 1, A, s01, x1, il::io, y0_par.Edit());
          },
          [&]() {
            dot_rec(parallelism - 1, A, s10, x0, il::io, y1);
          },
          [&]() {
            dot_rec(parallelism - 1, A, s11, x1, il::io, y1_par.Edit());
          });

      tbb::parallel_invoke(
          [&]() {
            for (il::int_t k = 0; k < n00; ++k) {
              y0[k] += y0_par[k];
            }
          },
          [&]() {
            for (il::int_t k = 0; k < n10; ++k) {
              y1[k] += y1_par[k];
            }
          });
    } else {
      dot_rec(parallelism - 1, A, s00, x0, il::io, y0);
      dot_rec(parallelism - 1, A, s01, x1, il::io, y0);
      dot_rec(parallelism - 1, A, s10, x0, il::io, y1);
      dot_rec(parallelism - 1, A, s11, x1, il::io, y1);
    }

    return;
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
il::Array<T> dot(const il::HMatrix<T>& A, const il::Array<T>& x) {
  IL_EXPECT_FAST(A.size(1) == x.size());

  il::Array<T> y{A.size(0), 0.0};
  il::ArrayEdit<T> y_edit = y.Edit();

  // For a nb of threads less than 4^parallelism
  il::int_t parallelism = 0;//4
  dot_rec(parallelism, A, A.root(), x.view(), il::io, y.Edit());

  return y;

}

}  // namespace il