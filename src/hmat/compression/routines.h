//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/algorithmArray.h>
#include <il/linearAlgebra/dense/blas/blas.h>
#include <il/linearAlgebra/dense/blas/dot.h>
#include <il/linearAlgebra/dense/factorization/LU.h>
#include <il/linearAlgebra/dense/factorization/Singular.h>

#include <hmat/arrayFunctor/MatrixGenerator.h>

namespace il {


inline double frobeniusNorm(const il::Array2D<double> &A) {
  double ans = 0.0;
  for (il::int_t i1 = 0; i1 < A.size(1); ++i1) {
    for (il::int_t i0 = 0; i0 < A.size(0); ++i0) {
      ans += A(i0, i1) * A(i0, i1);
    }
  }
  return ans;
}

inline double frobeniusNorm(const il::Array2D<std::complex<double>> &A) {
  double ans = 0.0;
  for (il::int_t i1 = 0; i1 < A.size(1); ++i1) {
    for (il::int_t i0 = 0; i0 < A.size(0); ++i0) {
      ans += il::ipow<2>(A(i0, i1).real()) + il::ipow<2>(A(i0, i1).imag());
    }
  }
  return ans;
}

template <il::int_t p, typename T>
void residual_row(const il::MatrixGenerator<T> &M,
                  const il::Array2D<T> &A, const il::Array2D<T> &B,
                  il::Range range0, il::Range range1, il::int_t i0, il::int_t r,
                  il::io_t, il::Array2DEdit<T> row) {
  IL_EXPECT_FAST(row.size(0) == p);
  IL_EXPECT_FAST(row.size(1) == (range1.end - range1.begin) * p);
  IL_EXPECT_FAST(range0.begin <= i0 && i0 < range0.end);

  const il::int_t n1 = range1.end - range1.begin;
  M.set(i0, range1.begin, il::io, row);
  if (r >= 1) {
    il::blas(
        -1.0,
        A.view(il::Range{(i0 - range0.begin) * p, (i0 - range0.begin + 1) * p},
               il::Range{0, r * p}),
        B.view(il::Range{0, r * p}, il::Range{0, n1 * p}), 1.0, il::io, row);
  }
};

template <il::int_t p, typename T>
void residual_column(const il::MatrixGenerator<T> &M,
                     const il::Array2D<T> &A, const il::Array2D<T> &B,
                     il::Range range0, il::Range range1, il::int_t i1,
                     il::int_t r, il::io_t, il::Array2DEdit<T> column) {
  IL_EXPECT_FAST(column.size(1) == p);
  IL_EXPECT_FAST(column.size(0) == (range0.end - range0.begin) * p);
  IL_EXPECT_FAST(range1.begin <= i1 && i1 < range1.end);

  const il::int_t n0 = range0.end - range0.begin;
  M.set(range0.begin, i1, il::io, column);
  if (r >= 1) {
    il::blas(
        -1.0, A.view(il::Range{0, n0 * p}, il::Range{0, r * p}),
        B.view(il::Range{0, r * p},
               il::Range{(i1 - range1.begin) * p, (i1 - range1.begin + 1) * p}),
        1.0, il::io, column);
  }
};

template <il::int_t p, typename T>
il::int_t find_largest_singular_value(const il::Array2D<T> &row,
                                      il::Range range1,
                                      const il::Array<il::int_t> &i1_used) {
  IL_EXPECT_FAST(row.size(0) == p);
  IL_EXPECT_FAST(row.size(1) == (range1.end - range1.begin) * p);

  const il::int_t n1 = range1.end - range1.begin;

  il::int_t i1_search = -1;
  double largest_singular_value = 0.0;

  for (il::int_t i1 = range1.begin; i1 < range1.end; ++i1) {
    bool already_searched = false;
    for (il::int_t k = 0; k < i1_used.size(); ++k) {
      if (i1 == i1_used[k]) {
        already_searched = true;
      }
    }
    if (!already_searched) {
      // To optimize: We don't need to compute the full list of singular
      // values. Only the smallest singular value needs to be computed
      il::StaticArray2D<T, p, p> matrix{};
      for (il::int_t j1 = 0; j1 < p; ++j1) {
        for (il::int_t j0 = 0; j0 < p; ++j0) {
          matrix(j0, j1) = row(j0, (i1 - range1.begin) * p + j1);
        }
      }

      il::Status status{};
      il::StaticArray<double, p> singular_values =
          il::singularValues(matrix, il::io, status);
      status.AbortOnError();
      il::sort(il::io, singular_values);

      if (singular_values[0] > largest_singular_value) {
        i1_search = i1;
        largest_singular_value = singular_values[0];
      }
    }
  }
  return i1_search;
}

template <il::int_t p, typename T>
il::StaticArray2D<T, p, p> lowRankSubmatrix(
    const il::MatrixGenerator<T> &M, const il::Array2D<T> &A,
    const il::Array2D<T> &B, il::int_t i0, il::int_t i1, il::int_t r) {
  il::StaticArray2D<T, p, p> matrix{0.0};
  if (r >= 1) {
    il::Array2DEdit<T> reference_matrix = matrix.Edit();
    il::blas(1.0, A.view(il::Range{i0 * p, (i0 + 1) * p}, il::Range{0, r * p}),
             B.view(il::Range{0, r * p}, il::Range{i1 * p, (i1 + 1) * p}), 0.0,
             il::io, reference_matrix);
  }
  return matrix;
};


template <il::int_t p, typename T>
il::int_t searchI0(const il::Array2D<T> &A, il::Range range0,
                   il::Range range1, const il::Array<il::int_t> i0_used,
                   il::int_t i1, il::int_t rank) {
  const il::int_t n0 = range0.end - range0.begin;
  const il::int_t n1 = range1.end - range0.begin;

  il::int_t i0_search = -1;
  double largest_singular_value = 0.0;
  for (il::int_t i0 = range0.begin; i0 < range0.end; ++i0) {
    bool already_searched = false;
    for (il::int_t k = 0; k < i0_used.size(); ++k) {
      if (i0 == i0_used[k]) {
        already_searched = true;
      }
    }
    if (!already_searched) {
      // To optimize: We don't need to compute the full list of singular
      // values. Only the smallest singular value needs to be computed
      il::StaticArray2D<T, p, p> matrix{};
      for (il::int_t b1 = 0; b1 < p; ++b1) {
        for (il::int_t b0 = 0; b0 < p; ++b0) {
          matrix(b0, b1) = A((i0 - range0.begin) * p + b0, rank * p + b1);
        }
      }

      il::Status status{};
      il::StaticArray<double, p> singular_values =
          il::singularValues(matrix, il::io, status);
      status.AbortOnError();

      il::sort(il::io, singular_values);

      if (singular_values[0] > largest_singular_value) {
        i0_search = i0;
        largest_singular_value = singular_values[0];
      }
    }
  }
  return i0_search;
}

}  // namespace il
