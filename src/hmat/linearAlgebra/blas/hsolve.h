#pragma once

#ifdef IL_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#define IL_CBLAS_INT MKL_INT
#define IL_CBLAS_LAYOUT CBLAS_LAYOUT
#elif IL_OPENBLAS
#include <OpenBLAS/cblas.h>
#include <OpenBLAS/lapacke.h>
#define IL_CBLAS_INT int
#define IL_CBLAS_LAYOUT CBLAS_ORDER
#endif

//#include <Matrix.h>

#include <il/Array2DView.h>
#include <il/ArrayView.h>
#include <il/linearAlgebra/Matrix.h>
#include <il/linearAlgebra/dense/blas/solve.h>

#include <hmat/arrayFunctor/MatrixGenerator.h>
#include <hmat/linearAlgebra/blas/hblas.h>
#include <hmat/linearAlgebra/blas/hsolve.h>
#include <hmat/hmatrix/HMatrix.h>

#ifdef IL_PARALLEL
#include <tbb/tbb.h>
#endif

namespace il {

template <typename T>
void solve(const il::HMatrix<T>& lu, il::MatrixType type, il::io_t,
           il::ArrayEdit<T> xy);

template <typename T>
void solve(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
           il::ArrayEdit<T> x);

template <typename T>
void solveLower(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::ArrayEdit<T> x);

template <typename T>
void solveLower(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::Array2DEdit<T> A);

template <typename T>
void solveLower(double epsilon, const il::HMatrix<T>& lu, il::spot_t slu,
                il::spot_t s, il::io_t, il::HMatrix<T>& A);

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::ArrayEdit<T> x);

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::Array2DEdit<T> A);

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t slu, il::spot_t s,
                il::io_t, il::HMatrix<T>& A);

template <typename T>
void solveUpperTranspose(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                         il::Array2DEdit<T> A);

template <typename T>
void solveUpperTranspose(const il::HMatrix<T>& lu, il::spot_t slu, il::spot_t s,
                         il::io_t, il::HMatrix<T>& A);

template <typename T>
void solveUpperRight(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                     il::Array2DEdit<T> A);

template <typename T>
void solveUpperRight(double epsilon, const il::HMatrix<T>& lu, il::spot_t slu,
                     il::spot_t s, il::io_t, il::HMatrix<T>& A);

template <typename T>
void solve(const il::HMatrix<T>& lu, il::io_t, il::ArrayEdit<T> xy) {
  il::solve(lu, lu.root(), il::io, xy);
}

template <typename T>
void solve(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
           il::ArrayEdit<T> x) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == lu.size(1, s));
  IL_EXPECT_MEDIUM(lu.size(1, s) == x.size());

  solveLower(lu, s, il::io, x);
  solveUpper(lu, s, il::io, x);
}

template <typename T>
void solveLower(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::ArrayEdit<T> x) {
  il::Array2DEdit<T> x_as_matrix{x.Data(), x.size(), 1, x.size(), 0, 0};
  solveLower(lu, s, il::io, x_as_matrix);
}

template <typename T>
void solveLower(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::Array2DEdit<T> A) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == lu.size(1, s));
  IL_EXPECT_MEDIUM(lu.size(0, s) == A.size(0));

  if (lu.isFullLu(s)) {
    il::ArrayView<int> pivot = lu.asFullLuPivot(s);
    il::Array2DView<T> full = lu.asFullLu(s);
    il::solve(pivot, full, il::MatrixType::LowerUnit, il::io, A);
  } else if (lu.isHierarchical(s)) {
    const il::spot_t s00 = lu.child(s, 0, 0);
    const il::spot_t s10 = lu.child(s, 1, 0);
    const il::spot_t s11 = lu.child(s, 1, 1);
    const il::int_t n0 = lu.size(0, s00);
    const il::int_t n1 = lu.size(0, s11);
    solveLower(lu, s00, il::io,
               A.Edit(il::Range{0, n0}, il::Range{0, A.size(1)}));
    il::blas(T{-1.0}, lu, s10, il::MatrixType::Regular,
             A.view(il::Range{0, n0}, il::Range{0, A.size(1)}), T{1.0}, il::io,
             A.Edit(il::Range{n0, n0 + n1}, il::Range{0, A.size(1)}));
    solveLower(lu, s11, il::io,
               A.Edit(il::Range{n0, n0 + n1}, il::Range{0, A.size(1)}));
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void solveLower(double epsilon, const il::HMatrix<T>& lu, il::spot_t slu,
                il::spot_t s, il::io_t, il::HMatrix<T>& A) {
  IL_EXPECT_MEDIUM(lu.size(0, slu) == lu.size(1, slu));
  IL_EXPECT_MEDIUM(lu.size(1, slu) == A.size(0, s));

  if (lu.isFullLu(slu)) {
    il::ArrayView<int> pivot = lu.asFullLuPivot(slu);
    il::Array2DView<T> lower = lu.asFullLu(slu);
    if (A.isFullRank(s)) {
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::solve(pivot, lower, il::MatrixType::LowerUnit, il::io, full);
    } else if (A.isLowRank(s)) {
      il::Array2DEdit<T> full = A.AsLowRankA(s);
      il::solve(pivot, lower, il::MatrixType::LowerUnit, il::io, full);
    } else {
      IL_EXPECT_MEDIUM(A.isHierarchical(s));
      il::abort();
      // Not sure if this needs to be implemented
    }
  } else if (lu.isHierarchical(slu)) {
    const il::spot_t slu00 = lu.child(slu, 0, 0);
    const il::spot_t slu10 = lu.child(slu, 1, 0);
    const il::spot_t slu11 = lu.child(slu, 1, 1);
    const il::int_t n0 = lu.size(0, slu00);
    const il::int_t n1 = lu.size(0, slu10);
    if (A.isFullRank(s)) {
      // FIXME: Are we sure we really need this il::MatrixType::Regular
      // parameter?
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::Array2DEdit<T> full0 =
          full.Edit(il::Range{0, n0}, il::Range{0, full.size(1)});
      il::Array2DEdit<T> full1 =
          full.Edit(il::Range{n0, n0 + n1}, il::Range{0, full.size(1)});
      il::solveLower(lu, slu00, il::io, full0);
      il::blas(T{-1.0}, lu, slu10, il::MatrixType::Regular, full0, T{1.0},
               il::io, full1);
      il::solveLower(lu, slu11, il::io, full1);
    } else if (A.isLowRank(s)) {
      il::Array2DEdit<T> lowA = A.AsLowRankA(s);
      il::Array2DEdit<T> lowA0 =
          lowA.Edit(il::Range{0, n0}, il::Range{0, lowA.size(1)});
      il::Array2DEdit<T> lowA1 =
          lowA.Edit(il::Range{n0, n0 + n1}, il::Range{0, lowA.size(1)});
      il::solveLower(lu, slu00, il::io, lowA0);
      il::blas(T{-1.0}, lu, slu10, il::MatrixType::Regular, lowA0, T{1.0},
               il::io, lowA1);
      il::solveLower(lu, slu11, il::io, lowA1);
    } else {
      const il::spot_t s00 = A.child(s, 0, 0);
      const il::spot_t s01 = A.child(s, 0, 1);
      const il::spot_t s10 = A.child(s, 1, 0);
      const il::spot_t s11 = A.child(s, 1, 1);
#ifdef IL_PARALLEL
      tbb::parallel_invoke(
          [&] { il::solveLower(epsilon, lu, slu00, s00, il::io, A); },
          [&] { il::solveLower(epsilon, lu, slu00, s01, il::io, A); });
      tbb::parallel_invoke(
          [&] {
            il::blas(epsilon, T{-1.0}, lu, slu10, A, s00, T{1.0}, s10, il::io,
                     A);
          },
          [&] {
            il::blas(epsilon, T{-1.0}, lu, slu10, A, s01, T{1.0}, s11, il::io,
                     A);
          });
      tbb::parallel_invoke(
          [&] { il::solveLower(epsilon, lu, slu11, s10, il::io, A); },
          [&] { il::solveLower(epsilon, lu, slu11, s11, il::io, A); });
#else
      il::solveLower(epsilon, lu, slu00, s00, il::io, A);
      il::solveLower(epsilon, lu, slu00, s01, il::io, A);
      il::blas(epsilon, T{-1.0}, lu, slu10, A, s00, T{1.0}, s10, il::io, A);
      il::blas(epsilon, T{-1.0}, lu, slu10, A, s01, T{1.0}, s11, il::io, A);
      il::solveLower(epsilon, lu, slu11, s10, il::io, A);
      il::solveLower(epsilon, lu, slu11, s11, il::io, A);
#endif
    }
  } else {
    IL_UNREACHABLE;
    il::abort();
  }
}

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::ArrayEdit<T> x) {
  il::Array2DEdit<T> x_as_matrix{x.Data(), x.size(), 1, x.size(), 0, 0};
  solveUpper(lu, s, il::io, x_as_matrix);
}

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                il::Array2DEdit<T> A) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == lu.size(1, s));
  IL_EXPECT_MEDIUM(lu.size(1, s) == A.size(0));

  if (lu.isFullLu(s)) {
    il::Array2DView<T> full_lu = lu.asFullLu(s);
    il::solve(full_lu, il::MatrixType::UpperNonUnit, il::io, A);
  } else if (lu.isHierarchical(s)) {
    const il::spot_t s00 = lu.child(s, 0, 0);
    const il::spot_t s01 = lu.child(s, 0, 1);
    const il::spot_t s11 = lu.child(s, 1, 1);
    const il::int_t n0 = lu.size(0, s00);
    const il::int_t n1 = lu.size(0, s11);
    solveUpper(lu, s11, il::io,
               A.Edit(il::Range{n0, n0 + n1}, il::Range{0, A.size(1)}));
    il::blas(T{-1.0}, lu, s01, il::MatrixType::Regular,
             A.view(il::Range{n0, n0 + n1}, il::Range{0, A.size(1)}), T{1.0},
             il::io, A.Edit(il::Range{0, n0}, il::Range{0, A.size(1)}));
    solveUpper(lu, s00, il::io,
               A.Edit(il::Range{0, n0}, il::Range{0, A.size(1)}));
  } else {
    IL_UNREACHABLE;
    il::abort();
  }
}

template <typename T>
void solveUpper(const il::HMatrix<T>& lu, il::spot_t slu, il::spot_t s,
                il::io_t, il::HMatrix<T>& A) {
  IL_EXPECT_MEDIUM(lu.size(0, slu) == lu.size(1, slu));
  IL_EXPECT_MEDIUM(lu.size(1, slu) == A.size(0, s));

  if (lu.isFullRank(slu)) {
    il::ArrayView<int> pivot = lu.asFullLuPivot(slu);
    il::Array2DView<T> upper = lu.asFullRank(slu);
    if (A.isFullRank(s)) {
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::solve(upper, il::MatrixType::UpperNonUnit, il::io, full);
    } else if (A.isLowRank(s)) {
      il::Array2DEdit<T> full = A.AsLowRankA(s);
      il::solve(upper, il::MatrixType::UpperNonUnit, il::io, full);
    } else {
      IL_EXPECT_MEDIUM(A.isHierarchical(s));
      IL_UNREACHABLE;
      il::abort();
      // Not sure if this needs to be implemented
    }
  } else if (lu.isHierarchical(slu)) {
    const il::spot_t s00 = lu.child(slu, 0, 0);
    const il::spot_t s01 = lu.child(slu, 0, 1);
    const il::spot_t s11 = lu.child(slu, 1, 1);
    const il::int_t n0 = lu.size(0, s00);
    const il::int_t n1 = lu.size(0, s11);
    if (A.isFullRank(s)) {
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::Array2DEdit<T> full0 =
          full.Edit(il::Range{0, n0}, il::Range{0, full.size(1)});
      il::Array2DEdit<T> full1 =
          full.Edit(il::Range{n0, n0 + n1}, il::Range{0, full.size(1)});

      il::solveUpper(lu, s11, il::io, full1);
      il::blas(T{-1.0}, lu, s01, il::MatrixType::Regular, full1, T{1.0}, il::io,
               full0);
      il::solveUpper(lu, s00, il::io, full0);
    } else if (A.isLowRank(s)) {
      il::Array2DEdit<T> lowA = A.AsLowRankA(s);
      il::Array2DEdit<T> lowA0 =
          lowA.Edit(il::Range{0, n0}, il::Range{0, lowA.size(1)});
      il::Array2DEdit<T> lowA1 =
          lowA.Edit(il::Range{n0, n0 + n1}, il::Range{0, lowA.size(1)});
      il::solveUpper(lu, s11, il::io, lowA1);
      il::blas(T{-1.0}, lu, s01, il::MatrixType::Regular, lowA1, T{1.0}, il::io,
               lowA0);
      il::solveLower(lu, s00, il::io, lowA0);
    } else {
      IL_EXPECT_MEDIUM(A.isHierarchical(s));
      IL_UNREACHABLE;
      // Not sure if this needs to be implemented
      il::abort();
    }
  } else {
    IL_UNREACHABLE;
  }
}

template <typename T>
void solveUpperTranspose(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                         il::Array2DEdit<T> A) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == lu.size(1, s));
  IL_EXPECT_MEDIUM(lu.size(0, s) == A.size(0));

  if (lu.isFullLu(s)) {
    il::Array2DView<T> full_lu = lu.asFullLu(s);
    il::solve(full_lu, il::MatrixType::UpperNonUnit, il::Dot::Transpose, il::io,
              A);
  } else if (lu.isHierarchical(s)) {
    const il::spot_t s00 = lu.child(s, 0, 0);
    const il::spot_t s01 = lu.child(s, 0, 1);
    const il::spot_t s11 = lu.child(s, 1, 1);
    const il::int_t n0 = lu.size(0, s00);
    const il::int_t n1 = lu.size(0, s11);
    il::Array2DEdit<T> A0 = A.Edit(il::Range{0, n0}, il::Range{0, A.size(1)});
    il::Array2DEdit<T> A1 =
        A.Edit(il::Range{n0, n0 + n1}, il::Range{0, A.size(1)});
    il::solveUpperTranspose(lu, s00, il::io, A0);
    il::blas(T{-1.0}, lu, s01, il::Dot::Transpose, A0, T{1.0}, il::io, A1);
    il::solveUpperTranspose(lu, s11, il::io, A1);
  } else {
    IL_UNREACHABLE;
    il::abort();
  }
}

template <typename T>
void solveUpperTranspose(const il::HMatrix<T>& lu, il::spot_t slu, il::spot_t s,
                         il::io_t, il::HMatrix<T>& A) {
  il::abort();
}

template <typename T>
void solveUpperRight(const il::HMatrix<T>& lu, il::spot_t s, il::io_t,
                     il::Array2DEdit<T> A) {
  IL_EXPECT_MEDIUM(lu.size(0, s) == lu.size(1, s));
  IL_EXPECT_MEDIUM(lu.size(0, s) == A.size(1));

  if (lu.isFullLu(s)) {
    // Validated
    il::Array2DView<T> full_lu = lu.asFullLu(s);
    il::solveRight(full_lu, il::MatrixType::UpperNonUnit, il::io, A);
  } else if (lu.isHierarchical(s)) {
    // Validated
    const il::spot_t s00 = lu.child(s, 0, 0);
    const il::spot_t s01 = lu.child(s, 0, 1);
    const il::spot_t s11 = lu.child(s, 1, 1);
    const il::int_t n0 = lu.size(0, s00);
    const il::int_t n1 = lu.size(0, s11);
    il::Array2DEdit<T> A0 = A.Edit(il::Range{0, A.size(0)}, il::Range{0, n0});
    il::Array2DEdit<T> A1 =
        A.Edit(il::Range{0, A.size(0)}, il::Range{n0, n0 + n1});
    solveUpperRight(lu, s00, il::io, A0);
    il::blas(T{-1.0}, A0, lu, s01, T{1.0}, il::io, A1);
    solveUpperRight(lu, s11, il::io, A1);
  } else {
    IL_UNREACHABLE;
    il::abort();
  }
}

template <typename T>
void solveUpperRight(double epsilon, const il::HMatrix<T>& lu, il::spot_t slu,
                     il::spot_t s, il::io_t, il::HMatrix<T>& A) {
  IL_EXPECT_MEDIUM(lu.size(0, slu) == lu.size(1, slu));
  IL_EXPECT_MEDIUM(lu.size(0, slu) == A.size(1, s));

  if (lu.isFullLu(slu)) {
    il::Array2DView<T> upper = lu.asFullLu(slu);
    if (A.isFullRank(s)) {
      // validated
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::solveRight(upper, il::MatrixType::UpperNonUnit, il::io, full);
    } else if (A.isLowRank(s)) {
      // Validated
      il::Array2DEdit<T> ab = A.AsLowRankB(s);
      il::solve(upper, il::MatrixType::UpperNonUnit, il::Dot::Transpose, il::io,
                ab);
    } else {
      IL_EXPECT_MEDIUM(A.isHierarchical(s));
      IL_UNREACHABLE;
      il::abort();
      // Not sure if this needs to be implemented
    }
  } else if (lu.isHierarchical(slu)) {
    const il::spot_t slu00 = lu.child(slu, 0, 0);
    const il::spot_t slu01 = lu.child(slu, 0, 1);
    const il::spot_t slu11 = lu.child(slu, 1, 1);
    const il::int_t n0 = lu.size(0, slu00);
    const il::int_t n1 = lu.size(0, slu11);
    if (A.isFullRank(s)) {
      // Validated
      il::Array2DEdit<T> full = A.AsFullRank(s);
      il::Array2DEdit<T> full0 =
          full.Edit(il::Range{0, full.size(0)}, il::Range{0, n0});
      il::Array2DEdit<T> full1 =
          full.Edit(il::Range{0, full.size(0)}, il::Range{n0, n0 + n1});

      il::solveUpperRight(lu, slu00, il::io, full0);
      il::blas(T{-1.0}, full0, lu, slu01, T{1.0}, il::io, full1);
      il::solveUpperRight(lu, slu11, il::io, full1);
    } else if (A.isLowRank(s)) {
      // Validated
      il::Array2DEdit<T> lowb = A.AsLowRankB(s);
      il::Array2DEdit<T> lowb0 =
          lowb.Edit(il::Range{0, n0}, il::Range{0, lowb.size(1)});
      il::Array2DEdit<T> lowb1 =
          lowb.Edit(il::Range{n0, n0 + n1}, il::Range{0, lowb.size(1)});
      il::solveUpperTranspose(lu, slu00, il::io, lowb0);
      il::blas(T{-1.0}, lu, slu01, il::Dot::Transpose, lowb0, T{1.0}, il::io,
               lowb1);
      il::solveUpperTranspose(lu, slu11, il::io, lowb1);
    } else {
      IL_EXPECT_MEDIUM(A.isHierarchical(s));
      // Validated
      const il::spot_t s00 = A.child(s, 0, 0);
      const il::spot_t s01 = A.child(s, 0, 1);
      const il::spot_t s10 = A.child(s, 1, 0);
      const il::spot_t s11 = A.child(s, 1, 1);
#ifdef IL_PARALLEL
      tbb::parallel_invoke(
          [&] { il::solveUpperRight(epsilon, lu, slu00, s00, il::io, A); },
          [&] { il::solveUpperRight(epsilon, lu, slu00, s10, il::io, A); });
      tbb::parallel_invoke(
          [&] {
            il::blas(epsilon, T{-1.0}, A, s00, lu, slu01, T{1.0}, s01, il::io,
                     A);
          },
          [&] {
            il::blas(epsilon, T{-1.0}, A, s10, lu, slu01, T{1.0}, s11, il::io,
                     A);
          });
      tbb::parallel_invoke(
          [&] { il::solveUpperRight(epsilon, lu, slu11, s01, il::io, A); },
          [&] { il::solveUpperRight(epsilon, lu, slu11, s11, il::io, A); });
#else
      il::solveUpperRight(epsilon, lu, slu00, s00, il::io, A);
      il::solveUpperRight(epsilon, lu, slu00, s10, il::io, A);
      il::blas(epsilon, T{-1.0}, A, s00, lu, slu01, T{1.0}, s01, il::io, A);
      il::blas(epsilon, T{-1.0}, A, s10, lu, slu01, T{1.0}, s11, il::io, A);
      il::solveUpperRight(epsilon, lu, slu11, s01, il::io, A);
      il::solveUpperRight(epsilon, lu, slu11, s11, il::io, A);
#endif
    }
  } else {
    IL_UNREACHABLE;
    il::abort();
  }
}
}  // namespace il