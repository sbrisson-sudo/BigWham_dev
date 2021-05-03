
#include <gtest/gtest.h>

#ifdef IL_BLAS
//#include <il/algorithmArray.h>
//#include <il/linearAlgebra/dense/blas/dot.h>
//#include <il/linearAlgebra/dense/blas/blas_static.h>
//#include <il/linearAlgebra/dense/factorization/Singular.h>
//
//TEST(Singular, test1) {
//  il::StaticArray2D<double, 3, 3> A{
//      il::value, {{3.0, 0.0, 0.0}, {0.0, 5.0, 0.0}, {0.0, 0.0, 7.0}}};
//  const double theta = 1.0;
//  A = il::dot(A, rotation1(theta));
//  A = il::dot(rotation2(2 * theta), A);
//  A = il::dot(A, rotation3(3 * theta));
//
//  il::Status status{};
//  il::StaticArray<double, 3> singular_value =
//      il::singularValues(A, il::io, status);
//  status.AbortOnError();
//  il::sort(il::io, singular_value);
//
//  const double epsilon = 1.0e-15;
//
//  ASSERT_TRUE(il::abs(singular_value[0] - 3.0) <= epsilon &&
//              il::abs(singular_value[1] - 5.0) <= epsilon &&
//              il::abs(singular_value[2] - 7.0) <= epsilon);
//}

#endif  // IL_BLAS
