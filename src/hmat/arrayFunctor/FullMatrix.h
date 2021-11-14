//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#include <il/Array2D.h>
#include <il/math.h>

#include <Hmat-lib/arrayFunctor/MatrixGenerator.h>

namespace il {

template <typename T>
class FullMatrix : public MatrixGenerator<T> {
 private:
  il::Array2D<T> A_;

 public:
  FullMatrix(il::Array2D<T> A);
  il::int_t size(il::int_t d) const override;
  il::int_t blockSize() const override;
  il::int_t sizeAsBlocks(il::int_t d) const override;
  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<T> M) const override;
};

template <typename T>
FullMatrix<T>::FullMatrix(il::Array2D<T> A) : A_{std::move(A)} {};

template <typename T>
il::int_t FullMatrix<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return A_.size(d);
};

template <typename T>
il::int_t FullMatrix<T>::blockSize() const {
  return 1;
}

template <typename T>
il::int_t FullMatrix<T>::sizeAsBlocks(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return A_.size(d);
}

template <typename T>
void FullMatrix<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                        il::Array2DEdit<T> M) const {
  IL_EXPECT_MEDIUM(b0 + M.size(0) <= size(0));
  IL_EXPECT_MEDIUM(b1 + M.size(1) <= size(1));

  for (il::int_t i1 = 0; i1 < M.size(1); ++i1) {
    for (il::int_t i0 = 0; i0 < M.size(0); ++i0) {
      M(i0, i1) = A_(b0 + i0, b1 + i1);
    }
  }
}

}  // namespace il
