#pragma once

#include <iostream>

#include <il/linearAlgebra/dense/factorization/luDecomposition.h>

#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/linearAlgebra/blas/hsolve.h>

#ifdef IL_PARALLEL
#include <tbb/tbb.h>
#endif

namespace il {

template <typename T>
inline void luDecomposition(double epsilon, il::io_t, il::HMatrix<T>& H);

template <typename T>
 inline void luDecomposition(double epsilon, il::spot_t s, il::io_t, il::HMatrix<T>& H);

template <typename T>
inline void luDecomposition(double epsilon, il::io_t, il::HMatrix<T>& H) {
  luDecomposition(epsilon, H.root(), il::io, H);
}

template <typename T>
inline void luDecomposition(double epsilon, il::spot_t s, il::io_t,
                     il::HMatrix<T>& H) {
  if (H.isFullRank(s)) {
    H.ConvertToFullLu(s);
    il::ArrayEdit<int> pivot = H.AsFullLuPivot(s);
    il::Array2DEdit<T> lu = H.AsFullLu(s);
    il::luDecomposition(il::io, pivot, lu);
  } else if (H.isHierarchical(s)) {
    const il::spot_t s00 = H.child(s, 0, 0);
    const il::spot_t s01 = H.child(s, 0, 1);
    const il::spot_t s10 = H.child(s, 1, 0);
    const il::spot_t s11 = H.child(s, 1, 1);
    il::Timer timer00{};
    timer00.Start();
    il::luDecomposition(epsilon, s00, il::io, H);
    timer00.Stop();
#ifdef IL_PARALLEL
    tbb::parallel_invoke(
        [&] { il::solveLower(epsilon, H, s00, s01, il::io, H); },
        [&] { il::solveUpperRight(epsilon, H, s00, s10, il::io, H); });
#else
    il::Timer timer01{};
    timer01.Start();
    il::solveLower(epsilon, H, s00, s01, il::io, H);
    timer01.Stop();
    il::Timer timer10{};
    timer10.Start();
    il::solveUpperRight(epsilon, H, s00, s10, il::io, H);
    timer10.Stop();
#endif
    il::Timer timer11{};
    timer11.Start();
    il::blas(epsilon, T{-1.0}, H, s10, H, s01, T{1.0}, s11, il::io, H);
    il::luDecomposition(epsilon, s11, il::io, H);
    timer11.Stop();
//    if (s.index == 0) {
//      std::cout << "Time00: " << timer00.time() << "s" << std::endl;
#ifdef IL_PARALLEL
#else
//      std::cout << "Time01: " << timer01.time() << "s" << std::endl;
//      std::cout << "Time10: " << timer10.time() << "s" << std::endl;
#endif
//      std::cout << "Time11: " << timer11.time() << "s" << std::endl;
//    }
  }
}

}  // namespace il
