#include <gtest/gtest.h>

//#include <hmat/hmatrix/HMatrixUtils.h>
//
//TEST(toArray2D, test0) {
//  const il::int_t n = 2;
//
//  il::HMatrix<double> H{};
//  il::spot_t s = H.root();
//  H.SetHierarchical(s);
//  H.SetFullRank(H.child(s, 0, 0), 2, 2);
//  H.SetFullRank(H.child(s, 1, 1), 2, 2);
//  H.SetLowRank(H.child(s, 1, 0), 2, 2, 1);
//  H.SetLowRank(H.child(s, 0, 1), 2, 2, 1);
//
//  {
//    il::Array2DEdit<double> E = H.AsFullRank(H.child(s, 0, 0));
//    E(0, 0) = 1.0;
//    E(1, 0) = 2.0;
//    E(0, 1) = 3.0;
//    E(1, 1) = 4.0;
//  }
//
//  {
//    il::Array2DEdit<double> E = H.AsFullRank(H.child(s, 1, 1));
//    E(0, 0) = 5.0;
//    E(1, 0) = 6.0;
//    E(0, 1) = 7.0;
//    E(1, 1) = 8.0;
//  }
//
//  {
//    il::Array2DEdit<double> EA = H.AsLowRankA(H.child(s, 1, 0));
//    il::Array2DEdit<double> EB = H.AsLowRankB(H.child(s, 1, 0));
//    EA(0, 0) = 1.0;
//    EA(1, 0) = 0.5;
//    EB(0, 0) = 2.0;
//    EB(1, 0) = 4.0;
//  }
//
//  {
//    il::Array2DEdit<double> EA = H.AsLowRankA(H.child(s, 0, 1));
//    il::Array2DEdit<double> EB = H.AsLowRankB(H.child(s, 0, 1));
//    EA(0, 0) = 0.25;
//    EA(1, 0) = 0.5;
//    EB(0, 0) = 8.0;
//    EB(1, 0) = 16.0;
//  }
//
//  il::Array2D<double> A = il::toArray2D(H);
//
//  ASSERT_TRUE(A.size(0) == 4 && A.size(1) == 4 && A(0, 0) == 1.0 &&
//              A(1, 0) == 2.0 && A(2, 0) == 2.0 && A(3, 0) == 1.0 &&
//              A(0, 1) == 3.0 && A(1, 1) == 4.0 && A(2, 1) == 4.0 &&
//              A(3, 1) == 2.0 && A(0, 2) == 2.0 && A(1, 2) == 4.0 &&
//              A(2, 2) == 5.0 && A(3, 2) == 6.0 && A(0, 3) == 4.0 &&
//              A(1, 3) == 8.0 && A(2, 3) == 7.0 && A(3, 3) == 8.0);
//}
