
#pragma once

#include <tbb/tbb.h>

#include <il/Array.h>
#include <il/linearAlgebra.h>

#include <hmat/hmatrix/HMatrix.h>

namespace il {

template <typename T>
void dot_rec(il::int_t parallelism, const bie::HMatrix<T>& A, il::spot_t s,
             il::ArrayView<T> x, il::io_t, il::ArrayEdit<T> y) {
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
          [&]() { dot_rec(parallelism - 1, A, s00, x0, il::io, y0); },
          [&]() {
            dot_rec(parallelism - 1, A, s01, x1, il::io, y0_par.Edit());
          },
          [&]() { dot_rec(parallelism - 1, A, s10, x0, il::io, y1); },
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
il::Array<T> dot(const bie::HMatrix<T>& A, const il::Array<T>& x) {
  IL_EXPECT_FAST(A.size(1) == x.size());

  il::Array<T> y{A.size(0), 0.0};
  // il::ArrayEdit<T> y_edit = y.Edit();

  // For a nb of threads less than 4^parallelism
  il::int_t parallelism = 4;  // 4
  dot_rec(parallelism, A, A.root(), x.view(), il::io, y.Edit());

  return y;
}



// full-rank dot product
static void dotfullrank(const bie::HMatrix<double>& A_,
                        const il::Array2D<il::int_t>& fullRank_pattern,
                        ArrayView<double> x_, il::io_t, il::Array<double>& y_) {
  IL_EXPECT_FAST(A_.size(1) == x_.size());
  IL_EXPECT_FAST(A_.size(0) == y_.size());

  /// serial loop
  for (il::int_t i = 0; i < fullRank_pattern.size(1); i++) {
    il::spot_t s(fullRank_pattern(0, i));
    il::int_t i0 = fullRank_pattern(1, i) - 1;
    il::int_t j0 = fullRank_pattern(2, i) - 1;

    il::Array2DView<double> a = A_.asFullRank(s);
    il::int_t i1 = i0 + a.size(0);
    il::int_t j1 = j0 + a.size(1);

    il::ArrayView<double> xs = x_.view(il::Range{j0, j1});
    il::ArrayEdit<double> ys = y_.Edit(il::Range{i0, i1});

    il::blas(1.0, a, xs, 1.0, il::io, ys);
  }
}

// low-rank block dot product
static void dotlowrank(const bie::HMatrix<double>& A_,
                       const il::Array2D<il::int_t>& lowRank_pattern,
                       ArrayView<double> x_, il::io_t, il::Array<double>& y_) {
  IL_EXPECT_FAST(A_.size(1) == x_.size());
  IL_EXPECT_FAST(A_.size(0) == y_.size());

  /// serial loop
  for (il::int_t i = 0; i < lowRank_pattern.size(1); i++) {
    il::spot_t s(lowRank_pattern(0, i));

    il::int_t i0 = lowRank_pattern(1, i) - 1;
    il::int_t j0 = lowRank_pattern(2, i) - 1;
    il::Array2DView<double> a = A_.asLowRankA(s);
    il::Array2DView<double> b = A_.asLowRankB(s);
    il::int_t i1 = i0 + a.size(0);
    il::int_t j1 = j0 + b.size(0);
    il::ArrayView<double> xs = x_.view(il::Range{j0, j1});
    il::ArrayEdit<double> ys = y_.Edit(il::Range{i0, i1});
    const il::int_t r = a.size(1);
    il::Array<double> tmp{r, 0.0};
    il::blas(1.0, b, il::Dot::Transpose, xs, 0.0, il::io, tmp.Edit());
    il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
  }
}

// dot product with pattern using tbb parallel_invoke
inline il::Array<double> dotwithpattern(
    const bie::HMatrix<double>& A_, const il::Array2D<il::int_t>& FR_pattern,
    const il::Array2D<il::int_t>& LR_pattern,const Array<double>& x_) {

  IL_EXPECT_FAST(A_.size(1) == x_.size());

  il::Array<double> y_FR{A_.size(0), 0.0};

  il::Array<double> y_LR{A_.size(0), 0.0};
   tbb::parallel_invoke([&] { dotfullrank(A_, FR_pattern, x_.view(), il::io, y_FR); },
                       [&] { dotlowrank(A_, LR_pattern, x_.view(), il::io, y_LR); });
   IL_EXPECT_FAST(y_FR.size() == y_LR.size());

  for (il::int_t i = 0; i < y_FR.size(); i++) {
    y_FR[i] += y_LR[i];
  }

  return y_FR;
}


//  serial dot product with the stored pattern
static il::Array<double> dotwithpattern_serial(
    const bie::HMatrix<double>& A_, const il::Array2D<il::int_t>& h_pattern,
    const il::Array<double>& x_) {
  IL_EXPECT_FAST(A_.size(1) == x_.size());

  il::Array<double> y_{A_.size(0), 0.0};

  /// serial loop
  for (il::int_t i = 0; i < h_pattern.size(1); i++) {
    il::spot_t s(h_pattern(0, i));
    il::int_t i0 = h_pattern(1, i) - 1;
    il::int_t j0 = h_pattern(2, i) - 1;

    if (A_.isFullRank(s)) {
      il::Array2DView<double> a = A_.asFullRank(s);
      il::int_t i1 = i0 + a.size(0);
      il::int_t j1 = j0 + a.size(1);

      il::ArrayView<double> xs = x_.view(il::Range{j0, j1});
      il::ArrayEdit<double> ys = y_.Edit(il::Range{i0, i1});

      il::blas(1.0, a, xs, 1.0, il::io, ys);

    } else if (A_.isLowRank(s)) {
      il::Array2DView<double> a = A_.asLowRankA(s);
      il::Array2DView<double> b = A_.asLowRankB(s);
      il::int_t i1 = i0 + a.size(0);
      il::int_t j1 = j0 + b.size(0);
      il::ArrayView<double> xs = x_.view(il::Range{j0, j1});
      il::ArrayEdit<double> ys = y_.Edit(il::Range{i0, i1});
      const il::int_t r = a.size(1);
      il::Array<double> tmp{r, 0.0};
      il::blas(1.0, b, il::Dot::Transpose, xs, 0.0, il::io, tmp.Edit());
      il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
    }
  }
  return y_;

}

}  // namespace il