#pragma once

#include <il/blas.h>

#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/hmatrix/LowRank.h>
#include <Hmat-lib/linearAlgebra/blas/hblas.h>
#include <Hmat-lib/linearAlgebra/factorization/lowRankApproximation.h>

namespace il {

template <typename T>
void blas(double epsilon, T alpha, const il::HMatrix<T>& A, il::spot_t sa,
          const il::HMatrix<T>& B, il::spot_t sb, T beta, il::spot_t sc,
          il::io_t, il::HMatrix<T>& C);

template <typename T>
void blas(T alpha, const il::HMatrix<T>& A, il::spot_t s, il::Array2DView<T> B,
          T beta, il::io_t, il::Array2DEdit<T> C);

template <typename T>
void blas(T alpha, const il::HMatrix<T>& A, il::spot_t s, il::Dot op,
          il::Array2DView<T> B, T beta, il::io_t, il::Array2DEdit<T> C);

template <typename T>
void blas(T alpha, il::Array2DView<T> A, const il::HMatrix<T>& B, il::spot_t s,
          T beta, il::io_t, il::Array2DEdit<T> C);

// Adds the Low Rank matrix A.B^T to the Hierachical matrix C
template <typename T>
void blasLowRank(double epsilon, T alpha, il::Array2DView<T> A,
                 il::Array2DView<T> B, T beta, il::spot_t s, il::io_t,
                 il::HMatrix<T>& C);

template <typename T>
void blas_rec(double alpha, const il::HMatrix<double>& A, il::spot_t s,
              il::MatrixType type, il::ArrayView<double> x, double beta,
              il::io_t, il::ArrayEdit<double> y);

template <typename T>
void blas(double alpha, const il::HMatrix<double>& lu, il::spot_t s,
          il::MatrixType type, il::ArrayView<double> x, double beta, il::io_t,
          il::ArrayEdit<double> y);

template <typename T>
void blas_rec(double alpha, const il::HMatrix<double>& A, il::spot_t s,
              il::MatrixType type, il::Array2DView<double> B, double beta,
              il::io_t, il::Array2DEdit<double> C);

template <typename T>
void blas(double alpha, const il::HMatrix<double>& lu, il::spot_t s,
          il::MatrixType type, il::Array2DView<double> A, double beta, il::io_t,
          il::Array2DEdit<double> B);

template <typename T>
void blas(double epsilon, T alpha, const il::HMatrix<T>& A, il::spot_t sa,
          const il::HMatrix<T>& B, il::spot_t sb, T beta, il::spot_t sc,
          il::io_t, il::HMatrix<T>& C) {
  IL_EXPECT_MEDIUM(A.size(0, sa) == C.size(0, sc));
  IL_EXPECT_MEDIUM(B.size(1, sb) == C.size(1, sc));
  IL_EXPECT_MEDIUM(A.size(1, sa) == B.size(0, sb));

  if (C.isFullRank(sc)) {
    il::Array2DEdit<T> c = C.AsFullRank(sc);
    if (A.isFullRank(sa) && B.isFullRank(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::Array2DView<T> b = B.asFullRank(sb);
      il::blas(alpha, a, b, beta, il::io, c);
    } else if (A.isFullRank(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp{a.size(0), ba.size(1)};
      il::blas(T{1.0}, a, ba, T{0.0}, il::io, tmp.Edit());
      il::blas(alpha, tmp.view(), bb, il::Dot::Transpose, beta, il::io, c);
    } else if (A.isLowRank(sa) && B.isFullRank(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2DView<T> b = A.asFullRank(sb);
      il::Array2D<T> tmp{ab.size(1), b.size(1)};
      il::blas(T{1.0}, ab, il::Dot::Transpose, b, T{0.0}, il::io, tmp.Edit());
      il::blas(alpha, aa, tmp.view(), beta, il::io, c);
    } else if (A.isLowRank(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp0{ab.size(1), ba.size(1)};
      il::blas(T{1.0}, ab, il::Dot::Transpose, ba, T{0.0}, il::io, tmp0.Edit());
      il::Array2D<T> tmp1{aa.size(0), ba.size(1)};
      il::blas(T{1.0}, aa, tmp0.view(), T{0.0}, il::io, tmp1.Edit());
      il::blas(alpha, tmp1.view(), bb, il::Dot::Transpose, beta, il::io, c);
    } else if (A.isFullRank(sa) && B.isHierarchical(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::blas(alpha, a, B, sb, beta, il::io, c);
    } else if (A.isHierarchical(sa) && B.isFullRank(sb)) {
      il::Array2DView<T> b = B.asFullRank(sb);
      il::blas(alpha, A, sa, b, beta, il::io, c);
    } else if (A.isLowRank(sa) && B.isHierarchical(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2D<T> tmp{B.size(1, sb), ab.size(1)};
      il::blas(T{1.0}, B, sb, il::Dot::Transpose, ab, T{0.0}, il::io,
               tmp.Edit());
      il::blas(alpha, aa, tmp.view(), il::Dot::Transpose, beta, il::io, c);
    } else if (A.isHierarchical(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp{A.size(0, sa), ba.size(1)};
      il::blas(T{1.0}, A, sa, ba, T{0.0}, il::io, tmp.Edit());
      il::blas(alpha, tmp.view(), bb, il::Dot::Transpose, beta, il::io, c);
    } else if (A.isHierarchical(sa) && B.isHierarchical(sb)) {
      il::LowRank<T> lrb = il::lowRank(epsilon, B, sb);
      il::Array2D<T> tmp{A.size(0, sa), lrb.A.size(1)};
      il::blas(T{1.0}, A, sa, lrb.A.view(), T{0.0}, il::io, tmp.Edit());
      il::blas(alpha, tmp.view(), lrb.B.view(), il::Dot::Transpose, beta,
               il::io, c);
    }
  } else if (C.isHierarchical(sc) || C.isLowRank(sc)) {
    if (A.isLowRank(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp0{ab.size(1), ba.size(1)};
      il::blas(T{1.0}, ab, il::Dot::Transpose, ba, T{0.0}, il::io, tmp0.Edit());
      il::Array2D<T> tmp1{aa.size(0), ba.size(1)};
      il::blas(T{1.0}, aa, tmp0.view(), T{0.0}, il::io, tmp1.Edit());
      il::blasLowRank(epsilon, alpha, tmp1.view(), bb, beta, sc, il::io, C);
    } else if (A.isFullRank(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp{a.size(0), ba.size(1)};
      il::blas(T{1.0}, a, ba, T{0.0}, il::io, tmp.Edit());
      il::blasLowRank(epsilon, alpha, tmp.view(), bb, beta, sc, il::io, C);
    } else if (A.isLowRank(sa) && B.isFullRank(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2DView<T> b = B.asFullRank(sb);
      il::Array2D<T> tmp{b.size(1), ab.size(1)};
      il::blas(T{1.0}, b, il::Dot::Transpose, ab, T{0.0}, il::io, tmp.Edit());
      il::blasLowRank(epsilon, alpha, aa, tmp.view(), beta, sc, il::io, C);
    } else if (A.isHierarchical(sa) && B.isLowRank(sb)) {
      il::Array2DView<T> ba = B.asLowRankA(sb);
      il::Array2DView<T> bb = B.asLowRankB(sb);
      il::Array2D<T> tmp{A.size(0, sa), ba.size(1)};
      il::blas(T{1.0}, A, sa, ba, T{0.0}, il::io, tmp.Edit());
      il::blasLowRank(epsilon, alpha, tmp.view(), bb, beta, sc, il::io, C);
    } else if (A.isLowRank(sa) && B.isHierarchical(sb)) {
      il::Array2DView<T> aa = A.asLowRankA(sa);
      il::Array2DView<T> ab = A.asLowRankB(sa);
      il::Array2D<T> tmp{B.size(1, sb), ab.size(1)};
      il::blas(T{1.0}, B, sb, il::Dot::Transpose, ab, T{0.0}, il::io,
               tmp.Edit());
      il::blasLowRank(epsilon, alpha, aa, tmp.view(), beta, sc, il::io, C);
    } else if (A.isFullRank(sa) && B.isFullRank(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::Array2DView<T> b = A.asFullRank(sb);
      il::Array2D<T> tmp0{a.size(0), b.size(1)};
      il::blas(T{1.0}, a, b, T{0.0}, il::io, tmp0.Edit());
      il::Array2D<T> tmp1{tmp0.size(1), tmp0.size(1), T{0.0}};
      for (il::int_t i = 0; i < tmp1.size(0); ++i) {
        tmp1(i, i) = T{1.0};
      }
      il::blasLowRank(epsilon, alpha, tmp0.view(), tmp1.view(), beta, sc,
                      il::io, C);
    } else if (A.isFullRank(sa) && B.isHierarchical(sb)) {
      il::Array2DView<T> a = A.asFullRank(sa);
      il::LowRank<T> lrb = il::lowRank(epsilon, B, sb);
      il::Array2D<T> tmp{a.size(0), lrb.A.size(1)};
      il::blas(T{1.0}, a, lrb.A.view(), T{0.0}, il::io, tmp.Edit());
      il::blasLowRank(epsilon, alpha, tmp.view(), lrb.B.view(), beta, sc,
                      il::io, C);
    } else if (A.isHierarchical(sa) && B.isFullRank(sb)) {
      il::LowRank<T> lra = il::lowRank(epsilon, A, sa);
      il::Array2DView<T> b = B.asFullRank(sb);
      il::Array2D<T> tmp{b.size(1), lra.B.size(1)};
      il::blas(T{1.0}, b, il::Dot::Transpose, lra.B.view(), T{0.0}, il::io,
               tmp.Edit());
      il::blasLowRank(epsilon, alpha, lra.A.view(), tmp.view(), beta, sc,
                      il::io, C);
    } else if (A.isHierarchical(sa) && B.isHierarchical(sb)) {
      if (C.isHierarchical(sc)) {
        il::blas(epsilon, alpha, A, A.child(sa, 0, 0), B, B.child(sb, 0, 0),
                 beta, C.child(sc, 0, 0), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 0, 1), B, B.child(sb, 1, 0),
                 T{1.0}, C.child(sc, 0, 0), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 0, 0), B, B.child(sb, 0, 1),
                 beta, C.child(sc, 0, 1), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 0, 1), B, B.child(sb, 1, 1),
                 T{1.0}, C.child(sc, 0, 1), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 1, 0), B, B.child(sb, 0, 0),
                 beta, C.child(sc, 1, 0), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 1, 1), B, B.child(sb, 1, 0),
                 T{1.0}, C.child(sc, 1, 0), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 1, 0), B, B.child(sb, 0, 1),
                 beta, C.child(sc, 1, 1), il::io, C);
        il::blas(epsilon, alpha, A, A.child(sa, 1, 1), B, B.child(sb, 1, 1),
                 T{1.0}, C.child(sc, 1, 1), il::io, C);
      } else {
        IL_EXPECT_FAST(C.isLowRank(sc));
        il::LowRank<T> lrb = il::lowRank(epsilon, B, sb);
        il::Array2D<T> tmp{A.size(0, sa), lrb.A.size(1)};
        il::blas(T{1.0}, A, sa, lrb.A.view(), T{0.0}, il::io, tmp.Edit());
        il::blasLowRank(epsilon, alpha, tmp.view(), lrb.B.view(), beta, sc,
                        il::io, C);
      }
    }
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas(T alpha, const il::HMatrix<T>& A, il::spot_t s, il::Array2DView<T> B,
          T beta, il::io_t, il::Array2DEdit<T> C) {
  IL_EXPECT_FAST(A.size(0, s) == C.size(0));
  IL_EXPECT_FAST(A.size(1, s) == B.size(0));
  IL_EXPECT_FAST(B.size(1) == C.size(1));

  if (A.isFullRank(s)) {
    il::Array2DView<T> a = A.asFullRank(s);
    il::blas(alpha, a, B, beta, il::io, C);
  } else if (A.isLowRank(s)) {
    il::Array2DView<T> aa = A.asLowRankA(s);
    il::Array2DView<T> ab = A.asLowRankB(s);
    const il::int_t r = aa.size(1);
    il::Array2D<T> tmp{r, B.size(1)};
    il::blas(T{1.0}, ab, il::Dot::Transpose, B, T{0.0}, il::io, tmp.Edit());
    il::blas(alpha, aa, tmp.view(), beta, il::io, C);
  } else if (A.isHierarchical(s)) {
    const il::spot_t s00 = A.child(s, 0, 0);
    const il::spot_t s10 = A.child(s, 1, 0);
    const il::spot_t s01 = A.child(s, 0, 1);
    const il::spot_t s11 = A.child(s, 1, 1);
    const il::int_t n00 = A.size(0, s00);
    const il::int_t n10 = A.size(0, s10);
    const il::int_t n01 = A.size(1, s00);
    const il::int_t n11 = A.size(1, s01);
    il::Array2DView<T> B0 = B.view(il::Range{0, n01}, il::Range{0, B.size(1)});
    il::Array2DView<T> B1 =
        B.view(il::Range{n01, n01 + n11}, il::Range{0, B.size(1)});
    il::Array2DEdit<T> C0 = C.Edit(il::Range{0, n00}, il::Range{0, C.size(1)});
    il::Array2DEdit<T> C1 =
        C.Edit(il::Range{n00, n00 + n10}, il::Range{0, C.size(1)});
    il::blas(alpha, A, s00, B0, beta, il::io, C0);
    il::blas(alpha, A, s01, B1, T{1.0}, il::io, C0);
    il::blas(alpha, A, s10, B0, beta, il::io, C1);
    il::blas(alpha, A, s11, B1, T{1.0}, il::io, C1);
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas(T alpha, const il::HMatrix<T>& A, il::spot_t s, il::Dot op,
          il::Array2DView<T> B, T beta, il::io_t, il::Array2DEdit<T> C) {
  IL_EXPECT_FAST(op == il::Dot::Transpose);
  IL_EXPECT_FAST(A.size(1, s) == C.size(0));
  IL_EXPECT_FAST(A.size(0, s) == B.size(0));
  IL_EXPECT_FAST(B.size(1) == C.size(1));

  if (A.isFullRank(s)) {
    il::Array2DView<T> a = A.asFullRank(s);
    il::blas(alpha, a, il::Dot::Transpose, B, beta, il::io, C);
  } else if (A.isLowRank(s)) {
    il::Array2DView<T> aa = A.asLowRankA(s);
    il::Array2DView<T> ab = A.asLowRankB(s);
    const il::int_t r = aa.size(1);
    il::Array2D<T> tmp{r, B.size(1)};
    il::blas(T{1.0}, aa, il::Dot::Transpose, B, T{0.0}, il::io, tmp.Edit());
    il::blas(alpha, ab, tmp.view(), beta, il::io, C);
  } else if (A.isHierarchical(s)) {
    const il::spot_t s00 = A.child(s, 0, 0);
    const il::spot_t s10 = A.child(s, 1, 0);
    const il::spot_t s01 = A.child(s, 0, 1);
    const il::spot_t s11 = A.child(s, 1, 1);
    const il::int_t n00 = A.size(0, s00);
    const il::int_t n10 = A.size(0, s10);
    const il::int_t n01 = A.size(1, s00);
    const il::int_t n11 = A.size(1, s01);
    il::Array2DView<T> B0 = B.view(il::Range{0, n00}, il::Range{0, B.size(1)});
    il::Array2DView<T> B1 =
        B.view(il::Range{n00, n00 + n10}, il::Range{0, B.size(1)});
    il::Array2DEdit<T> C0 = C.Edit(il::Range{0, n01}, il::Range{0, C.size(1)});
    il::Array2DEdit<T> C1 =
        C.Edit(il::Range{n01, n01 + n11}, il::Range{0, C.size(1)});
    il::blas(alpha, A, s00, il::Dot::Transpose, B0, beta, il::io, C0);
    il::blas(alpha, A, s10, il::Dot::Transpose, B1, T{1.0}, il::io, C0);
    il::blas(alpha, A, s01, il::Dot::Transpose, B0, beta, il::io, C1);
    il::blas(alpha, A, s11, il::Dot::Transpose, B1, T{1.0}, il::io, C1);
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas(T alpha, il::Array2DView<T> A, const il::HMatrix<T>& B, il::spot_t s,
          T beta, il::io_t, il::Array2DEdit<T> C) {
  IL_EXPECT_FAST(A.size(0) == C.size(0));
  IL_EXPECT_FAST(A.size(1) == B.size(0, s));
  IL_EXPECT_FAST(B.size(1, s) == C.size(1));

  if (B.isFullRank(s)) {
    il::Array2DView<T> b = B.asFullRank(s);
    il::blas(alpha, A, b, beta, il::io, C);
  } else if (B.isLowRank(s)) {
    il::Array2DView<T> ba = B.asLowRankA(s);
    il::Array2DView<T> bb = B.asLowRankB(s);
    const il::int_t r = ba.size(1);
    il::Array2D<T> tmp{A.size(0), r};
    il::blas(T{1.0}, A, ba, T{0.0}, il::io, tmp.Edit());
    il::blas(alpha, tmp.view(), bb, il::Dot::Transpose, beta, il::io, C);
  } else if (B.isHierarchical(s)) {
    const il::spot_t s00 = B.child(s, 0, 0);
    const il::spot_t s10 = B.child(s, 1, 0);
    const il::spot_t s01 = B.child(s, 0, 1);
    const il::spot_t s11 = B.child(s, 1, 1);
    const il::int_t n00 = B.size(0, s00);
    const il::int_t n10 = B.size(0, s10);
    const il::int_t n01 = B.size(1, s00);
    const il::int_t n11 = B.size(1, s01);
    il::Array2DView<T> A0 = A.view(il::Range{0, A.size(0)}, il::Range{0, n00});
    il::Array2DView<T> A1 =
        A.view(il::Range{0, A.size(0)}, il::Range{n00, n00 + n10});
    il::Array2DEdit<T> C0 = C.Edit(il::Range{0, C.size(0)}, il::Range{0, n01});
    il::Array2DEdit<T> C1 =
        C.Edit(il::Range{0, C.size(0)}, il::Range{n01, n01 + n11});
    il::blas(alpha, A0, B, s00, beta, il::io, C0);
    il::blas(alpha, A1, B, s10, beta, il::io, C0);
    il::blas(alpha, A0, B, s01, beta, il::io, C1);
    il::blas(alpha, A1, B, s11, beta, il::io, C1);
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blasLowRank(double epsilon, T alpha, il::Array2DView<T> A,
                 il::Array2DView<T> B, T beta, il::spot_t s, il::io_t,
                 il::HMatrix<T>& C) {
  IL_EXPECT_FAST(A.size(1) == B.size(1));
  IL_EXPECT_FAST(A.size(0) == C.size(0, s));
  IL_EXPECT_FAST(B.size(0) == C.size(1, s));

  if (C.isFullRank(s)) {
    il::Array2DEdit<T> c = C.AsFullRank(s);
    il::blas(alpha, A, B, il::Dot::Transpose, beta, il::io, c);
  } else if (C.isLowRank(s)) {
    const il::int_t r0 = C.rankOfLowRank(s);
    const il::int_t r1 = A.size(1);
    if (epsilon == T{0.0}) {
      C.UpdateRank(s, r0 + r1);
      il::Array2DEdit<T> ca = C.AsLowRankA(s);
      il::Array2DEdit<T> cb = C.AsLowRankB(s);
      il::Array2DEdit<T> cb_old =
          cb.Edit(il::Range{0, cb.size(0)}, il::Range{0, r0});
      for (il::int_t i1 = 0; i1 < cb_old.size(1); ++i1) {
        for (il::int_t i0 = 0; i0 < cb_old.size(0); ++i0) {
          cb_old(i0, i1) *= beta;
        }
      }
      il::Array2DEdit<T> ca_new =
          ca.Edit(il::Range{0, ca.size(0)}, il::Range{r0, r0 + r1});
      il::Array2DEdit<T> cb_new =
          cb.Edit(il::Range{0, cb.size(0)}, il::Range{r0, r0 + r1});
      il::copy(A, il::io, ca_new);
      il::copy(B, il::io, cb_new);
      for (il::int_t i1 = 0; i1 < cb_new.size(1); ++i1) {
        for (il::int_t i0 = 0; i0 < cb_new.size(0); ++i0) {
          cb_new(i0, i1) *= alpha;
        }
      }
    } else {
      il::LowRank<T> ab = il::lowRankAddition(epsilon, alpha, A, B, beta,
                                              C.asLowRankA(s), C.asLowRankB(s));
      C.UpdateRank(s, ab.A.size(1));
      il::Array2DEdit<T> ca = C.AsLowRankA(s);
      il::Array2DEdit<T> cb = C.AsLowRankB(s);
      il::copy(ab.A.view(), il::io, ca);
      il::copy(ab.B.view(), il::io, cb);
    }
  } else if (C.isHierarchical(s)) {
    const il::spot_t s00 = C.child(s, 0, 0);
    const il::spot_t s10 = C.child(s, 1, 0);
    const il::spot_t s01 = C.child(s, 0, 1);
    const il::spot_t s11 = C.child(s, 1, 1);
    const il::int_t n00 = C.size(0, s00);
    const il::int_t n10 = C.size(0, s10);
    const il::int_t n01 = C.size(1, s00);
    const il::int_t n11 = C.size(1, s01);
    const il::int_t r = A.size(1);
    il::blasLowRank(epsilon, alpha, A.view(il::Range{0, n00}, il::Range{0, r}),
                    B.view(il::Range{0, n01}, il::Range{0, r}), beta, s00,
                    il::io, C);
    il::blasLowRank(epsilon, alpha, A.view(il::Range{0, n00}, il::Range{0, r}),
                    B.view(il::Range{n01, n01 + n11}, il::Range{0, r}), beta,
                    s01, il::io, C);
    il::blasLowRank(
        epsilon, alpha, A.view(il::Range{n00, n00 + n10}, il::Range{0, r}),
        B.view(il::Range{0, n01}, il::Range{0, r}), beta, s10, il::io, C);
    il::blasLowRank(epsilon, alpha,
                    A.view(il::Range{n00, n00 + n10}, il::Range{0, r}),
                    B.view(il::Range{n01, n01 + n11}, il::Range{0, r}), beta,
                    s11, il::io, C);
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas_rec(T alpha, const il::HMatrix<T>& A, il::spot_t s,
              il::MatrixType type, il::ArrayView<T> x, T beta, il::io_t,
              il::ArrayEdit<T> y) {
  IL_EXPECT_FAST(A.size(0, s) == y.size());
  IL_EXPECT_FAST(A.size(1, s) == x.size());

  if (A.isFullRank(s)) {
    il::Array2DView<T> a = A.asFullRank(s);
    if (type == il::MatrixType::Regular) {
      il::blas(alpha, a, x, beta, il::io, y);
    } else if (type == il::MatrixType::LowerUnit) {
      il::abort();
      IL_UNREACHABLE;
    } else {
      il::abort();
      IL_UNREACHABLE;
    }
    return;
  } else if (A.isLowRank(s)) {
    IL_EXPECT_MEDIUM(type == il::MatrixType::Regular);
    il::Array2DView<T> a = A.asLowRankA(s);
    il::Array2DView<T> b = A.asLowRankB(s);
    const il::int_t r = a.size(1);
    il::Array<T> tmp{r, T{0.0}};
    il::blas(T{1.0}, b, il::Dot::Transpose, x, T{0.0}, il::io, tmp.Edit());
    il::blas(alpha, a, tmp.view(), beta, il::io, y);
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
    if (type == il::MatrixType::LowerUnit) {
      il::blas_rec(alpha, A, s00, il::MatrixType::LowerUnit, x0, beta, il::io,
                   y0);
      il::blas_rec(alpha, A, s10, il::MatrixType::Regular, x0, beta, il::io,
                   y1);
      il::blas_rec(alpha, A, s11, il::MatrixType::LowerUnit, x1, beta, il::io,
                   y1);
    } else {
      IL_UNREACHABLE;
    }
    return;
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas(T alpha, const il::HMatrix<T>& lu, il::spot_t s, il::MatrixType type,
          il::ArrayView<T> x, T beta, il::io_t, il::ArrayEdit<T> y) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == y.size());
  IL_EXPECT_MEDIUM(lu.size(1, s) == x.size());

  il::blas_rec(alpha, lu, s, type, x, beta, il::io, y);
}

template <typename T>
void blas_rec(T alpha, const il::HMatrix<T>& A, il::spot_t s,
              il::MatrixType type, il::Array2DView<T> B, T beta, il::io_t,
              il::Array2DEdit<T> C) {
  IL_EXPECT_FAST(A.size(0, s) == C.size(0));
  IL_EXPECT_FAST(A.size(1, s) == B.size(0));

  if (A.isFullRank(s)) {
    il::Array2DView<T> a = A.asFullRank(s);
    if (type == il::MatrixType::Regular) {
      il::blas(alpha, a, B, beta, il::io, C);
    } else if (type == il::MatrixType::LowerUnit) {
      il::abort();
      IL_UNREACHABLE;
    } else {
      il::abort();
      IL_UNREACHABLE;
    }
    return;
  } else if (A.isLowRank(s)) {
    IL_EXPECT_MEDIUM(type == il::MatrixType::Regular);
    il::Array2DView<T> a = A.asLowRankA(s);
    il::Array2DView<T> b = A.asLowRankB(s);
    il::Array2D<T> tmp{b.size(1), B.size(1), T{0.0}};
    il::blas(T{1.0}, b, il::Dot::Transpose, B, T{0.0}, il::io, tmp.Edit());
    il::blas(alpha, a, tmp.view(), beta, il::io, C);
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
    il::Array2DView<T> B0 = B.view(il::Range{0, n01}, il::Range{0, B.size(1)});
    il::Array2DView<T> B1 =
        B.view(il::Range{n01, n01 + n11}, il::Range{0, B.size(1)});
    il::Array2DEdit<T> C0 = C.Edit(il::Range{0, n00}, il::Range{0, C.size(1)});
    il::Array2DEdit<T> C1 =
        C.Edit(il::Range{n00, n00 + n10}, il::Range{0, C.size(1)});
    if (type == il::MatrixType::LowerUnit) {
      // FIXME: Are we really going there?????
      // What the point of this il::MatrixType::LowerUnit anyway???
      // Anyway, I believe that the coefficients are wrong
      il::abort();
      il::blas_rec(alpha, A, s00, il::MatrixType::LowerUnit, B0, beta, il::io,
                   C0);
      il::blas_rec(alpha, A, s10, il::MatrixType::Regular, B0, beta, il::io,
                   C1);
      il::blas_rec(alpha, A, s11, il::MatrixType::LowerUnit, B1, beta, il::io,
                   C1);
    } else if (type == il::MatrixType::Regular) {
      il::blas_rec(alpha, A, s00, il::MatrixType::Regular, B0, beta, il::io,
                   C0);
      il::blas_rec(alpha, A, s01, il::MatrixType::Regular, B1, T{1.0}, il::io,
                   C0);
      il::blas_rec(alpha, A, s10, il::MatrixType::Regular, B0, beta, il::io,
                   C1);
      il::blas_rec(alpha, A, s11, il::MatrixType::Regular, B1, T{1.0}, il::io,
                   C1);
    }
    return;
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void blas(T alpha, const il::HMatrix<T>& lu, il::spot_t s, il::MatrixType type,
          il::Array2DView<T> A, T beta, il::io_t, il::Array2DEdit<T> B) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == B.size(0));
  IL_EXPECT_MEDIUM(lu.size(1, s) == A.size(0));

  il::blas_rec(alpha, lu, s, type, A, beta, il::io, B);
}

}  // namespace il