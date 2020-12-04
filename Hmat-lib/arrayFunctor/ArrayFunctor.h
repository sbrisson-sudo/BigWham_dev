#pragma once

#include <il/ArrayView.h>

#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/linearAlgebra/blas/hdot.h>

namespace il {

template <typename T>
class ArrayFunctor {
 public:
  virtual il::int_t size(il::int_t d) const = 0;
  virtual void operator()(il::ArrayView<T> x, il::io_t,
                          il::ArrayEdit<T> y) const = 0;
};

template <typename T>
class HMatrixGenerator : public ArrayFunctor<T> {
 private:
  il::HMatrix<T> h_;

 public:
  HMatrixGenerator(il::HMatrix<T> h);
  il::int_t size(il::int_t d) const override;
  void operator()(il::ArrayView<T> x, il::io_t,
                  il::ArrayEdit<T> y) const override;
};

template <typename T>
HMatrixGenerator<T>::HMatrixGenerator(il::HMatrix<T> h) : h_{std::move(h)} {}

template <typename T>
il::int_t HMatrixGenerator<T>::size(il::int_t d) const {
  return h_.size(d);
}

template <typename T>
void HMatrixGenerator<T>::operator()(il::ArrayView<T> x, il::io_t,
                                     il::ArrayEdit<T> y) const {
  il::dot_rec(h_, h_.root(), x, il::io, y);
}

}  // namespace il
