#pragma once

#include <il/Array2D.h>

namespace il {

template <typename T>
struct LowRank {
  il::Array2D<T> A;
  il::Array2D<T> B;
};

}  // namespace il