#include <gtest/gtest.h>

#include <il/Array.h>

#include <hmat/arrayFunctor/FullMatrix.h>
#include <hmat/arrayFunctor/GaussianMatrix.h>
#include <hmat/hmatrix/HMatrixType.h>
#include <hmat/hmatrix/LowRank.h>
#include <hmat/compression/adaptiveCrossApproximation.h>
//
//TEST(adaptiveCrossApproximation, test0) {
//  const il::int_t n = 4;
//  const double alpha = 1.0;
//  const il::GaussianMatrix<double> G{2 * n, il::Range{0, n},
//                                     il::Range{n, 2 * n}, alpha};
//
//  il::Tree<il::SubHMatrix, 4> tree{};
//  const il::spot_t s = tree.root();
//  tree.Set(s, il::SubHMatrix{il::Range{0, n}, il::Range{0, n},
//                             il::HMatrixType::LowRank});
//
//  const double epsilon = 1.0e-4;
//  const il::HMatrix<double> H = il::toHMatrix(G, tree, epsilon);
//
//  const il::Array2D<double> M0 = il::toArray2D(G);
//  const il::Array2D<double> M1 = il::toArray2D(H);
//  il::Array2D<double> diff{n, n};
//  for (il::int_t i1 = 0; i1 < n; ++i1) {
//    for (il::int_t i0 = 0; i0 < n; ++i0) {
//      diff(i0, i1) = M0(i0, i1) - M1(i0, i1);
//    }
//  }
//
//  ASSERT_TRUE(true);
//}

TEST(adaptiveCrossApproximation, test1) {

  /*
  Python code used to generate this test matrix of rank 2 :

  # 1. Compute A = U V^T
  U = np.arange(0., 1., 1/12).reshape((6,2))
  V = np.arange(1., 2., 1/12).reshape((6,2))
  A = U @ V.T  # This creates a 6x6 matrix of rank 2
  print("A =\n", A)
  print(f"Shape: {A.shape}, Rank: {np.linalg.matrix_rank(A)}")

  # 2. Write the corresponding C++ code string
  cpp_code = """il::Array2D<double> A{
      il::value,
  {{""" + \
  "},\n     {".join(
          [", ".join(f"{A[j,i]:.6f}" for j in range(A.shape[1])) for i in range(A.shape[0])]
      ) + "}}};"
  print("\nC++ code:\n")
  print(cpp_code)
  
  */

  il::Array2D<double> A{
      il::value,
    {{0.090278, 0.437500, 0.784722, 1.131944, 1.479167, 1.826389},
     {0.104167, 0.506944, 0.909722, 1.312500, 1.715278, 2.118056},
     {0.118056, 0.576389, 1.034722, 1.493056, 1.951389, 2.409722},
     {0.131944, 0.645833, 1.159722, 1.673611, 2.187500, 2.701389},
     {0.145833, 0.715278, 1.284722, 1.854167, 2.423611, 2.993056},
     {0.159722, 0.784722, 1.409722, 2.034722, 2.659722, 3.284722}}};
  const il::int_t n = 6;
  const il::FullMatrix<double> matgen{A};

  il::Range range0{0, n};
  il::Range range1{0, n};

  // Here we cerate a LowRank struct that stores the two Array2D we want to copy
  double epsilon = 1e-6;
  auto lra = bigwham::adaptiveCrossApproximation<1>(matgen, range0, range1, epsilon);

  // Compute matvec y=Ax  
  auto a = lra->A.view();
  auto b = lra->B.view();

  // std::cout << "rank = " << a.size(1) << "\n";

  il::Array<double> x{n, 0.};
  il::Array<double> y{n, 0.};
  auto x_edit = x.Edit();
  for (int i(0); i<n; i++) x_edit[i] = i; 

  il::Array<double> tmp{a.size(1), 0.};

  il::blas(1.0, b, il::Dot::None, x.view(), 0.0, il::io, tmp.Edit()); 
  il::blas(1.0, a, tmp.view(), 1.0, il::io, y.Edit());

  // std::cout << "y = [";
  // for (int i(0); i<n; i++) std::cout << y[i] << ", ";
  // std::cout << "]\n"; 

  // matvec result 
  std::vector<double> y_ref{2.11805556, 10.38194444, 18.64583333, 26.90972222, 35.17361111, 43.4375};
  double tol = 1e-4;

  for (size_t i = 0; i < n; ++i) {
    ASSERT_NEAR(y[i], y_ref[i], tol) << "Mismatch at index " << i;
  }

}
