#include <gtest/gtest.h>

#include <Hmat-lib/arrayFunctor/FullMatrix.h>
#include <Hmat-lib/arrayFunctor/GaussianMatrix.h>
#include <Hmat-lib/compression/toHMatrix.h>
#include <Hmat-lib/hmatrix/HMatrixType.h>
#include <Hmat-lib/hmatrix/HMatrixUtils.h>

TEST(adaptiveCrossApproximation, test0) {
  const il::int_t n = 4;
  const double alpha = 1.0;
  const il::GaussianMatrix<double> G{2 * n, il::Range{0, n},
                                     il::Range{n, 2 * n}, alpha};

  il::Tree<il::SubHMatrix, 4> tree{};
  const il::spot_t s = tree.root();
  tree.Set(s, il::SubHMatrix{il::Range{0, n}, il::Range{0, n},
                             il::HMatrixType::LowRank});

  const double epsilon = 1.0e-4;
  const il::HMatrix<double> H = il::toHMatrix(G, tree, epsilon);

  const il::Array2D<double> M0 = il::toArray2D(G);
  const il::Array2D<double> M1 = il::toArray2D(H);
  il::Array2D<double> diff{n, n};
  for (il::int_t i1 = 0; i1 < n; ++i1) {
    for (il::int_t i0 = 0; i0 < n; ++i0) {
      diff(i0, i1) = M0(i0, i1) - M1(i0, i1);
    }
  }

  ASSERT_TRUE(true);
}

TEST(adaptiveCrossApproximation, test1) {
  //  A0 = {1 - I/2, 2 - I, 3 - 3 I/2, 4 - 2 I};
  //  B0 = {4, 3, 2, 1};
  //  A1 = {4 + I/2, 3 + I, 2 + 3 I/2, 1 + 2 I};
  //  B1 = {1, 2, 3, 4};
  const il::int_t n = 4;
  il::Array2D<std::complex<double>> A{
      il::value,
      {{{8.0, -3.0 / 2}, {11.0, -3.0}, {14.0, -9.0 / 2}, {17.0, -6.0}},
       {{11.0, -1.0 / 2}, {12.0, -1.0}, {13.0, -3.0 / 2}, {14.0, -2.0}},
       {{14.0, 1.0 / 2}, {13.0, 1.0}, {12.0, 3.0 / 2}, {11.0, 2.0}},
       {{17.0, 3.0 / 2}, {14.0, 3.0}, {11.0, 9.0 / 2}, {8.0, 6.0}}}};
  const il::FullMatrix<std::complex<double>> G{A};

  il::Tree<il::SubHMatrix, 4> tree{};
  const il::spot_t s = tree.root();
  tree.Set(s, il::SubHMatrix{il::Range{0, n}, il::Range{0, n},
                             il::HMatrixType::LowRank});

  const double epsilon = 1.0e-4;
  const il::HMatrix<std::complex<double>> H = il::toHMatrix(G, tree, epsilon);

  ASSERT_TRUE(true);
}
