//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#ifndef BIGWHAM_GAUSSIANMATRIX_H
#define BIGWHAM_GAUSSIANMATRIX_H

#include <memory>

#include <il/math.h>

#include "hmat/arrayFunctor/matrix_generator.h"
#include "hmat/hierarchical_representation.h"

// This an example of a Matrix Generator needed to code for the implementation
// of a BIE kernel

namespace bigwham {

template <typename T> class GaussianMatrix : public bigwham::MatrixGenerator<T> {
private:
  il::int_t n_;
  il::Range range0_;
  il::Range range1_;
  double alpha_;

public:
  GaussianMatrix(il::int_t n, double alpha,
                 const std::shared_ptr<HRepresentation> &hr);
  GaussianMatrix(il::int_t n, il::Range range0, il::Range range1, double alpha,
                 const std::shared_ptr<HRepresentation> &hr);
  il::int_t size(il::int_t d) const override;
  il::int_t blockSize() const override;
  il::int_t sizeAsBlocks(il::int_t d) const override;
  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<double> M) const override;
};

template <typename T>
GaussianMatrix<T>::GaussianMatrix(il::int_t n, double alpha,
                                  const std::shared_ptr<HRepresentation> &hr) {
  IL_EXPECT_MEDIUM(n >= 0);

  this->hr_ = hr;
  n_ = n;
  range0_ = il::Range{0, n};
  range1_ = il::Range{0, n};
  alpha_ = alpha;
};

template <typename T>
GaussianMatrix<T>::GaussianMatrix(il::int_t n, il::Range range0,
                                  il::Range range1, double alpha,
                                  const std::shared_ptr<HRepresentation> &hr) {
  IL_EXPECT_MEDIUM(n >= 0);
  IL_EXPECT_MEDIUM(range0.begin >= 0 && range0.end <= n);
  IL_EXPECT_MEDIUM(range1.begin >= 0 && range1.end <= n);

  this->hr_ = hr;
  n_ = n;
  range0_ = range0;
  range1_ = range1;
  alpha_ = alpha;
};

template <typename T> il::int_t GaussianMatrix<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);
  switch (d) {
  case 0:
    return range0_.end - range0_.begin;
  case 1:
    return range1_.end - range1_.begin;
  default:
    IL_UNREACHABLE;
  }
  IL_UNREACHABLE;
  return -1;
};

template <typename T> il::int_t GaussianMatrix<T>::blockSize() const {
  return 1;
}

template <typename T>
il::int_t GaussianMatrix<T>::sizeAsBlocks(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return size(d);
}

template <typename T>
void GaussianMatrix<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                            il::Array2DEdit<double> M) const {
  IL_EXPECT_MEDIUM(b0 + M.size(0) <= size(0));
  IL_EXPECT_MEDIUM(b1 + M.size(1) <= size(1));

  const double beta = il::ipow<2>(1.0 / (n_ - 1));
  for (il::int_t i1 = 0; i1 < M.size(1); ++i1) {
    il::int_t j1 = range1_.begin + b1 + i1;
    for (il::int_t i0 = 0; i0 < M.size(0); ++i0) {
      il::int_t j0 = range0_.begin + b0 + i0;
      M(i0, i1) = std::exp(-alpha_ * il::ipow<2>(j0 - j1) * beta);
    }
  }
}

} // namespace bigwham

#endif