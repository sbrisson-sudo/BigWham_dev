//
// This file is part of BigWham
//
// Created by Brice Lecampion on 01.02.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

// this files contains Array Functor defining pre-conditioner and Matrix
// interface used for performing matrix free iterative solve

#ifndef BIGWHAM_HITERATIVESOLVER_H
#define BIGWHAM_HITERATIVESOLVER_H

#include <il/Array.h>
#include <il/Gmres.h>

#include <hmat/linearAlgebra/blas/hdot.h>

namespace bie {

/////////
template <typename T>
class DiagPreconditioner : public il::FunctorArray<T> {
  // diag preconditionner for Hmat gmres solve - note the diag array must have
  // been permuted here
 private:
  il::Array<T> diag_;

 public:
  explicit DiagPreconditioner(il::Array<T> diag) : diag_{std::move(diag)} {}

  il::int_t size(il::int_t d) const override {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);
    return diag_.size();
  }

  void operator()(il::ArrayView<T> x, il::io_t,
                  il::ArrayEdit<T> y) const override {
    const il::int_t n = diag_.size();
    IL_EXPECT_FAST(x.size() == n);
    IL_EXPECT_FAST(y.size() == n);

    for (il::int_t i = 0; i < n; ++i) {
      y[i] = diag_[i] * x[i];
    }
  }
};
/////////

/////////  here - must implement the proper permutation corresponding to Hlu
// template <typename T>
// class HLUPreconditioner : public il::FunctorArray<T> {
//  bie::HMatrix<T> Hlu_;
// public:
//  explicit HLUPreconditioner(bie::HMatrix<T>& Hlu) : Hlu_{std::move(Hlu)} {};
//
//  il::int_t size(il::int_t d) const override {
//    IL_EXPECT_MEDIUM(d == 0 || d == 1);
//    if (d == 0) {
//      return Hlu_.size(0);
//
//    } else if (d == 1) {
//      return Hlu_.size(1);
//    }
//  }
//
//  void operator()(il::ArrayView<T> x, il::io_t,
//                  il::ArrayEdit<T> y) const override {
//    const il::int_t rows = Hlu_.size(0);
//    const il::int_t columns = Hlu_.size(1);
//
//    IL_EXPECT_FAST(columns == rows);
//    IL_EXPECT_FAST(x.size() == columns);
//    IL_EXPECT_FAST(y.size() == rows);
//
//    for (il::int_t i = 0; i < rows; ++i) {
//      y[i] = x[i];
//    }
//    il::solve(Hlu_, il::io, y);
//  }
//};

/////////
template <typename T>
class Matrix : public il::FunctorArray<T> {
 private:
  bie::HMatrix<T> H_;
  il::Array2D<il::int_t>  fr_patt_;
  il::Array2D<il::int_t>  lr_patt_;

 public:
  explicit Matrix(bie::HMatrix<T>& H,il::Array2D<il::int_t> & fr_patt, il::Array2D<il::int_t> & lr_patt)
  : H_{H},fr_patt_{fr_patt},lr_patt_{lr_patt} {};

  il::int_t size(il::int_t d) const override {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);
    if (d == 0) {
      return H_.size(0);

    } else if (d == 1) {
      return H_.size(1);
    }
  }

  void operator()(il::ArrayView<T> x, il::io_t,
                  il::ArrayEdit<T> y) const override {
    const il::int_t rows = H_.size(0);
    const il::int_t columns = H_.size(1);

    IL_EXPECT_FAST(columns == rows);
    IL_EXPECT_FAST(x.size() == columns);
    IL_EXPECT_FAST(y.size() == rows);

    il::Array<double> z{columns};

    for (il::int_t i = 0; i < columns; i++) {
      z[i] = x[i];
    }

    z = il::dotwithpattern(H_,fr_patt_,lr_patt_,z);
    // Hmat product vector dot(H,x)   x must be an Array   NOT of ArrayView Type
    // hence the hack using an aux array
    // but we use memory...

    for (il::int_t i = 0; i < columns; i++) {
      y[i] = z[i];
    }
  }
};

}  // namespace bie

#endif  // BIGWHAM_HITERATIVESOLVER_H
