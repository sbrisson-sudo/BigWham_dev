#include <gtest/gtest.h>

#include <iostream>

#include <il/Tree.h>
#include <il/math.h>

#include <Hmat-lib/arrayFunctor/ArrayFunctor.h>
#include <Hmat-lib/cluster/cluster.h>
#include <Hmat-lib/compression/toHMatrix.h>
#include <Hmat-lib/hmatrix/HMatrix.h>
#include <Hmat-lib/hmatrix/HMatrixType.h>
#include <Hmat-lib/hmatrix/HMatrixUtils.h>
#include <Hmat-lib/linearAlgebra/blas/hdot.h>
#include <Hmat-lib/linearAlgebra/blas/hsolve.h>
#include <Hmat-lib/linearAlgebra/factorization/luDecomposition.h>
#include <Hmat-lib/arrayFunctor/GaussianMatrix.h>
#include <Hmat-lib/Matrix.h>      // Change this ....

TEST(lowRank, test0) {
  il::Array2D<double> A{il::value,
                        {{1.0, 0.0, 0.0, 0.0},
                         {0.0, 2.0, 0.0, 0.0},
                         {0.0, 0.0, 0.5, 0.0},
                         {0.0, 0.0, 0.0, 0.0}}};

  const double epsilon = 0.2;
  il::LowRank<double> ab = il::lowRank(epsilon, A.view());

  const int dummy = 0;
  (void)dummy;

  ASSERT_TRUE(true);
}

TEST(solve, test0) {
  il::HMatrix<double> H{};
  const il::spot_t s = H.root();
  H.SetHierarchical(s);

  const il::spot_t s00 = H.child(s, 0, 0);
  H.SetFullRank(s00, 2, 2);
  il::Array2DEdit<double> H00 = H.AsFullRank(s00);
  H00(0, 0) = 3.0;
  H00(1, 1) = 4.0;
  H00(0, 1) = 0.0;
  H00(1, 0) = 1.0;

  const il::spot_t s11 = H.child(s, 1, 1);
  H.SetFullRank(s11, 2, 2);
  il::Array2DEdit<double> H11 = H.AsFullRank(s11);
  H11(0, 0) = 5.0;
  H11(1, 1) = 6.0;
  H11(0, 1) = 0.0;
  H11(1, 0) = 0.0;

  const il::spot_t s01 = H.child(s, 0, 1);
  H.SetLowRank(s01, 2, 2, 1);
  il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
  il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
  H01A(0, 0) = 1.0;
  H01A(1, 0) = 1.0;
  H01B(0, 0) = 1.0;
  H01B(1, 0) = 1.0;

  const il::spot_t s10 = H.child(s, 1, 0);
  H.SetLowRank(s10, 2, 2, 1);
  il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
  il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
  H10A(0, 0) = 1.0;
  H10A(1, 0) = 1.0;
  H10B(0, 0) = 1.0;
  H10B(1, 0) = 1.0;

  const double epsilon_lu = 0.0;
  il::luDecomposition(epsilon_lu, il::io, H);

  il::Array<double> y = {il::value, {10.0, 16.0, 18.0, 27.0}};
  il::solve(H, il::io, y.Edit());

  const double eps = 1.0e-10;

  ASSERT_TRUE(il::abs(y[0] - 1.0) <= eps && il::abs(y[1] - 2.0) <= eps &&
              il::abs(y[2] - 3.0) <= eps && il::abs(y[3] - 4.0) <= eps);
}

TEST(solve, test1) {
  il::HMatrix<double> H{};
  const il::spot_t s = H.root();
  H.SetHierarchical(s);

  const il::spot_t s00 = H.child(s, 0, 0);
  H.SetFullRank(s00, 2, 2);
  il::Array2DEdit<double> H00 = H.AsFullRank(s00);
  H00(0, 0) = 1.0;
  H00(0, 1) = 4.0;
  H00(1, 0) = 3.0;
  H00(1, 1) = 0.0;

  const il::spot_t s11 = H.child(s, 1, 1);
  H.SetFullRank(s11, 2, 2);
  il::Array2DEdit<double> H11 = H.AsFullRank(s11);
  H11(0, 0) = 0.0;
  H11(0, 1) = 6.0;
  H11(1, 0) = 5.0;
  H11(1, 1) = 1.0;

  const il::spot_t s01 = H.child(s, 0, 1);
  H.SetLowRank(s01, 2, 2, 1);
  il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
  il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
  H01A(0, 0) = 1.0;
  H01A(1, 0) = 1.0;
  H01B(0, 0) = 1.0;
  H01B(1, 0) = 0.5;

  const il::spot_t s10 = H.child(s, 1, 0);
  H.SetLowRank(s10, 2, 2, 1);
  il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
  il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
  H10A(0, 0) = 0.5;
  H10A(1, 0) = 1.0;
  H10B(0, 0) = 1.0;
  H10B(1, 0) = 1.0;

  const double epsilon_lu = 0.0;
  il::luDecomposition(epsilon_lu, il::io, H);

  il::Array<double> y = {il::value, {14.0, 8.0, 25.5, 22.0}};
  il::solve(H, il::io, y.Edit());

  const double eps = 1.0e-14;

  ASSERT_TRUE(il::abs(y[0] - 1.0) <= eps && il::abs(y[1] - 2.0) <= eps &&
              il::abs(y[2] - 3.0) <= eps && il::abs(y[3] - 4.0) <= eps);
}

TEST(solve, test2) {
  il::HMatrix<double> H{};
  const il::spot_t s = H.root();
  H.SetHierarchical(s);

  {
    const il::spot_t s01 = H.child(s, 0, 1);
    H.SetLowRank(s01, 4, 4, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    for (il::int_t k = 0; k < 4; ++k) {
      H01A(k, 0) = 1.0;
      H01B(k, 0) = 1.0;
    }
  }

  {
    const il::spot_t s10 = H.child(s, 1, 0);
    H.SetLowRank(s10, 4, 4, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    for (il::int_t k = 0; k < 4; ++k) {
      H10A(k, 0) = 1.0;
      H10B(k, 0) = 1.0;
    }
  }

  {
    const il::spot_t s00outer = H.child(s, 0, 0);
    H.SetHierarchical(s00outer);

    const il::spot_t s00 = H.child(s00outer, 0, 0);
    H.SetFullRank(s00, 2, 2);
    il::Array2DEdit<double> H00 = H.AsFullRank(s00);
    H00(0, 0) = 10.0;
    H00(1, 1) = 11.0;
    H00(0, 1) = 1.0;
    H00(1, 0) = 1.0;

    const il::spot_t s11 = H.child(s00outer, 1, 1);
    H.SetFullRank(s11, 2, 2);
    il::Array2DEdit<double> H11 = H.AsFullRank(s11);
    H11(0, 0) = 12.0;
    H11(1, 1) = 13.0;
    H11(0, 1) = 1.0;
    H11(1, 0) = 1.0;

    const il::spot_t s01 = H.child(s00outer, 0, 1);
    H.SetLowRank(s01, 2, 2, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    H01A(0, 0) = 1.0;
    H01A(1, 0) = 1.0;
    H01B(0, 0) = 1.0;
    H01B(1, 0) = 1.0;

    const il::spot_t s10 = H.child(s00outer, 1, 0);
    H.SetLowRank(s10, 2, 2, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    H10A(0, 0) = 1.0;
    H10A(1, 0) = 1.0;
    H10B(0, 0) = 1.0;
    H10B(1, 0) = 1.0;
  }

  {
    const il::spot_t s11outer = H.child(s, 1, 1);
    H.SetHierarchical(s11outer);

    const il::spot_t s00 = H.child(s11outer, 0, 0);
    H.SetFullRank(s00, 2, 2);
    il::Array2DEdit<double> H00 = H.AsFullRank(s00);
    H00(0, 0) = 14.0;
    H00(1, 1) = 15.0;
    H00(0, 1) = 1.0;
    H00(1, 0) = 1.0;

    const il::spot_t s11 = H.child(s11outer, 1, 1);
    H.SetFullRank(s11, 2, 2);
    il::Array2DEdit<double> H11 = H.AsFullRank(s11);
    H11(0, 0) = 16.0;
    H11(1, 1) = 17.0;
    H11(0, 1) = 1.0;
    H11(1, 0) = 1.0;

    const il::spot_t s01 = H.child(s11outer, 0, 1);
    H.SetLowRank(s01, 2, 2, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    H01A(0, 0) = 1.0;
    H01A(1, 0) = 1.0;
    H01B(0, 0) = 1.0;
    H01B(1, 0) = 1.0;

    const il::spot_t s10 = H.child(s11outer, 1, 0);
    H.SetLowRank(s10, 2, 2, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    H10A(0, 0) = 1.0;
    H10A(1, 0) = 1.0;
    H10B(0, 0) = 1.0;
    H10B(1, 0) = 1.0;
  }

  il::Array<double> x_initial{H.size(1), 1.0};
  il::Array<double> y = il::dot(H, x_initial);

  const double epsilon_lu = 0.0;
  il::luDecomposition(epsilon_lu, il::io, H);

  il::solve(H, il::io, y.Edit());

  bool result = true;
  const double eps = 1.0e-14;
  for (il::int_t i = 0; i < y.size(); ++i) {
    if (il::abs(y[i] - 1.0) >= eps) {
      result = false;
    }
  }

  ASSERT_TRUE(result);
}

TEST(solve, test3) {
  il::HMatrix<double> H{};
  const il::spot_t s = H.root();
  H.SetHierarchical(s);

  {
    const il::spot_t s01 = H.child(s, 0, 1);
    H.SetLowRank(s01, 4, 4, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    for (il::int_t k = 0; k < 4; ++k) {
      H01A(k, 0) = 1.0 / (1 + 2 * k);
      H01B(k, 0) = 1.0 / (1 + k);
    }
  }

  {
    const il::spot_t s10 = H.child(s, 1, 0);
    H.SetLowRank(s10, 4, 4, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    for (il::int_t k = 0; k < 4; ++k) {
      H10A(k, 0) = 1.0 / (2 + k * k);
      H10B(k, 0) = 1.0 / (1 + k * k);
    }
  }

  {
    const il::spot_t s00outer = H.child(s, 0, 0);
    H.SetHierarchical(s00outer);

    const il::spot_t s00 = H.child(s00outer, 0, 0);
    H.SetFullRank(s00, 2, 2);
    il::Array2DEdit<double> H00 = H.AsFullRank(s00);
    H00(0, 0) = 10.0;
    H00(1, 1) = 11.0;
    H00(0, 1) = 1.0 / 3.14159;
    H00(1, 0) = 1.0 / 1.414;

    const il::spot_t s11 = H.child(s00outer, 1, 1);
    H.SetFullRank(s11, 2, 2);
    il::Array2DEdit<double> H11 = H.AsFullRank(s11);
    H11(0, 0) = 12.0;
    H11(1, 1) = 13.0;
    H11(0, 1) = 1.0 / 4.14159;
    H11(1, 0) = 1.0 / 2.414;

    const il::spot_t s01 = H.child(s00outer, 0, 1);
    H.SetLowRank(s01, 2, 2, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    H01A(0, 0) = 1.0 / 2;
    H01A(1, 0) = 1.0 / 7;
    H01B(0, 0) = 1.0 / 15;
    H01B(1, 0) = 1.0 / 12;

    const il::spot_t s10 = H.child(s00outer, 1, 0);
    H.SetLowRank(s10, 2, 2, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    H10A(0, 0) = 1.0 / 13;
    H10A(1, 0) = 1.0 / 14;
    H10B(0, 0) = 1.0 / 25;
    H10B(1, 0) = 1.0 / 31;
  }

  {
    const il::spot_t s11outer = H.child(s, 1, 1);
    H.SetHierarchical(s11outer);

    const il::spot_t s00 = H.child(s11outer, 0, 0);
    H.SetFullRank(s00, 2, 2);
    il::Array2DEdit<double> H00 = H.AsFullRank(s00);
    H00(0, 0) = 14.0;
    H00(1, 1) = 15.0;
    H00(0, 1) = 1.0 / 6.14159;
    H00(1, 0) = 1.0 / 4.414;

    const il::spot_t s11 = H.child(s11outer, 1, 1);
    H.SetFullRank(s11, 2, 2);
    il::Array2DEdit<double> H11 = H.AsFullRank(s11);
    H11(0, 0) = 16.0;
    H11(1, 1) = 17.0;
    H11(0, 1) = 1.0 / 7.14159;
    H11(1, 0) = 1.0 / 5.414;

    const il::spot_t s01 = H.child(s11outer, 0, 1);
    H.SetLowRank(s01, 2, 2, 1);
    il::Array2DEdit<double> H01A = H.AsLowRankA(s01);
    il::Array2DEdit<double> H01B = H.AsLowRankB(s01);
    H01A(0, 0) = 1.0 / 27;
    H01A(1, 0) = 1.0 / 23;
    H01B(0, 0) = 1.0 / 21;
    H01B(1, 0) = 1.0 / 42;

    const il::spot_t s10 = H.child(s11outer, 1, 0);
    H.SetLowRank(s10, 2, 2, 1);
    il::Array2DEdit<double> H10A = H.AsLowRankA(s10);
    il::Array2DEdit<double> H10B = H.AsLowRankB(s10);
    H10A(0, 0) = 1.0 / 35;
    H10A(1, 0) = 1.0 / 39;
    H10B(0, 0) = 1.0 / 41;
    H10B(1, 0) = 1.0 / 44;
  }

  il::Array<double> x_initial{H.size(1)};
  for (il::int_t i = 0; i < x_initial.size(); ++i) {
    x_initial[i] = static_cast<double>(1 + i);
  }
  il::Array<double> y = il::dot(H, x_initial);

  const double epsilon_lu = 0.0;
  il::luDecomposition(epsilon_lu, il::io, H);

  il::solve(H, il::io, y.Edit());

  bool result = true;
  const double eps = 1.0e-14;
  for (il::int_t i = 0; i < y.size(); ++i) {
    if (il::abs(y[i] - (1 + i)) >= eps) {
      result = false;
    }
  }

  ASSERT_TRUE(result);
}

TEST(solve, test4) {
  const il::int_t n = 32768;
  //  const il::int_t n = 8192;
  const il::int_t dim = 2;
  const il::int_t leaf_max_size = 256;

  const double radius = 1.0;
  il::Array2D<double> point{n, dim};
  for (il::int_t i = 0; i < n; ++i) {
    point(i, 0) = radius * std::cos((il::pi * (i + 0.5)) / n);
    point(i, 1) = radius * std::sin((il::pi * (i + 0.5)) / n);
  }
  const il::Cluster cluster = il::cluster(leaf_max_size, il::io, point);

  //////////////////////////////////////////////////////////////////////////////
  // Here we prepare the compression scheme of the matrix depending upon the
  // clustering and the relative diameter and distance in between the clusters
  //////////////////////////////////////////////////////////////////////////////
  //
  // And then, the clustering of the matrix, only from the geometry of the
  // points.
  // - Use eta = 0 for no compression
  // - Use eta = 1 for moderate compression
  // - Use eta = 10 for large compression
  // We have compression when Max(diam0, diam1) <= eta * distance
  // Use a large value for eta when you want everything to be Low-Rank when
  // you are outside the diagonal
//  const double eta = 10.0;
  const double eta = 0.1;
  const il::Tree<il::SubHMatrix, 4> hmatrix_tree =
      il::hmatrixTree(point, cluster.partition, eta);

  //////////////////////////////////////////////////////////////////////////////
  // We build the H-Matrix
  //////////////////////////////////////////////////////////////////////////////
  const double alpha = 10.0;
  // Put mu / n, and the condition number should be mu
  const double beta = 0.1 / n;
  const il::Matrix<double> M{point, alpha, beta};
  const double epsilon = 0.1;
  il::HMatrix<double> h = il::toHMatrix(M, hmatrix_tree, epsilon);

  // First, we compute the compression ratio
  std::cout << "Compression ratio: " << il::compressionRatio(h) << std::endl;

  il::Array<double> x{h.size(0), 1.0};
  il::Array<double> y = il::dot(h, x);

  //////////////////////////////////////////////////////////////////////////////
  // We convert it to a regular matrix to compute its condition number
  //////////////////////////////////////////////////////////////////////////////
  //  const il::Array2D<double> full_h = il::toArray2D(h);
  //  il::Status status{};
  il::Timer timer{};
  //  timer.Start();
  //  const il::LU<il::Array2D<double>> full_lu_h{full_h, il::io, status};
  //  timer.Stop();
  //  status.AbortOnError();
  //
  //  std::cout << "Time for full LU-decomposition: " << timer.time() <<
  //  std::endl;
  //
  //  const double norm_full_h = il::norm(full_h, il::Norm::L1);
  //  const double cn = full_lu_h.conditionNumber(il::Norm::L1, norm_full_h);
  //
  //  il::Array<double> y_full = full_lu_h.solve(y);
  //
  //  double relative_error_full = 0.0;
  //  for (il::int_t i = 0; i < y_full.size(); ++i) {
  //    const double re = il::abs(y_full[i] - 1.0);
  //    if (re > relative_error_full) {
  //      relative_error_full = re;
  //    }
  //  }

  //////////////////////////////////////////////////////////////////////////////
  // Let's play with it
  //////////////////////////////////////////////////////////////////////////////
  const double epsilon_lu = 1.0e-10;
  //  const double epsilon_lu = 0;
  timer.Reset();
  timer.Start();
  il::luDecomposition(epsilon_lu, il::io, h);
  timer.Stop();
  il::solve(h, il::io, y.Edit());

  std::cout << "Time for HLU-decomposition: " << timer.time() << std::endl;
  std::cout << "Compression ratio of HLU-decomposition: "
            << il::compressionRatio(h) << std::endl;

  double relative_error = 0.0;
  for (il::int_t i = 0; i < y.size(); ++i) {
    const double re = il::abs(y[i] - 1.0);
    if (re > relative_error) {
      relative_error = re;
    }
  }

  ASSERT_TRUE(relative_error <= 1.0e-10);
}

TEST(solve, test5) {
  const il::int_t n = 32768;
  const il::int_t dim = 2;
  const il::int_t leaf_max_size = 256;

  const double radius = 1.0;
  il::Array2D<double> point{n, dim};
  for (il::int_t i = 0; i < n; ++i) {
    point(i, 0) = radius * std::cos((il::pi * (i + 0.5)) / n);
    point(i, 1) = radius * std::sin((il::pi * (i + 0.5)) / n);
  }
  const il::Cluster cluster = il::cluster(leaf_max_size, il::io, point);

    const double eta = 0.1;
//  const double eta = 10.0;
  const il::Tree<il::SubHMatrix, 4> hmatrix_tree =
      il::hmatrixTree(point, cluster.partition, eta);

  const double alpha = 10.0;
  const double beta = 0.1 / n;
  const il::Matrix<std::complex<double>> M{point, alpha, beta};
  const double epsilon = 0.1;
  il::HMatrix<std::complex<double>> h = il::toHMatrix(M, hmatrix_tree, epsilon);

  // First, we compute the compression ratio
  std::cout << "Compression ratio: " << il::compressionRatio(h) << std::endl;

  il::Array<std::complex<double>> x{h.size(0), 1.0};
  il::Array<std::complex<double>> y = il::dot(h, x);

  //////////////////////////////////////////////////////////////////////////////
  // Let's play with it
  //////////////////////////////////////////////////////////////////////////////
  const double epsilon_lu = 1.0e-10;
  //  const double epsilon_lu = 0;
  il::Timer timer{};
  timer.Start();
  il::luDecomposition(epsilon_lu, il::io, h);
  timer.Stop();
  il::solve(h, il::io, y.Edit());

  std::cout << "Time for HLU-decomposition: " << timer.time() << std::endl;
  std::cout << "Compression ratio of HLU-decomposition: "
            << il::compressionRatio(h) << std::endl;

  double relative_error = 0.0;
  for (il::int_t i = 0; i < y.size(); ++i) {
    const double re = il::abs(y[i] - 1.0);
    if (re > relative_error) {
      relative_error = re;
    }
  }

  ASSERT_TRUE(relative_error <= 1.0e-10);
}
