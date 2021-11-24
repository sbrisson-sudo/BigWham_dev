#pragma once

#include <hmat/hmatrix/HMatrix.h>
#include <hmat/hmatrix/LowRank.h>
#include <complex>

namespace il {

bie::LowRank<double> lowRank(double epsilon, il::Array2DView<double> A);

bie::LowRank<std::complex<double>> lowRank(
    double epsilon, il::Array2DView<std::complex<double>> A);

template <typename T>
bie::LowRank<T> lowRank(double epsilon, const bie::HMatrix<T>& A, il::spot_t s);

bie::LowRank<double> lowRankAddition(double epsilon, double alpha,
                                    il::Array2DView<double> aa,
                                    il::Array2DView<double> ab, double beta,
                                    il::Array2DView<double> ba,
                                    il::Array2DView<double> bb);

bie::LowRank<std::complex<double>> lowRankAddition(
    double epsilon, std::complex<double> alpha,
    il::Array2DView<std::complex<double>> aa,
    il::Array2DView<std::complex<double>> ab, std::complex<double> beta,
    il::Array2DView<std::complex<double>> ba,
    il::Array2DView<std::complex<double>> bb);

template <typename T>
bie::LowRank<T> lowRank(double epsilon, const bie::HMatrix<T>& A, il::spot_t s) {
  if (A.isFullRank(s)) {
    return il::lowRank(epsilon, A.asFullRank(s));
  } else if (A.isLowRank(s)) {
    il::Array2DView<T> a = A.asLowRankA(s);
    il::Array2DView<T> b = A.asLowRankB(s);
    bie::LowRank<T> ans{};
    ans.A.Resize(a.size(0), a.size(1));
    il::copy(a, il::io, ans.A.Edit());
    ans.B.Resize(b.size(0), b.size(1));
    il::copy(b, il::io, ans.B.Edit());
    return ans;
  } else if (A.isHierarchical(s)) {
    const il::spot_t s00 = A.child(s, 0, 0);
    const il::spot_t s01 = A.child(s, 0, 1);
    const il::spot_t s10 = A.child(s, 1, 0);
    const il::spot_t s11 = A.child(s, 1, 1);
    bie::LowRank<T> lr00 = il::lowRank(epsilon, A, s00);
    bie::LowRank<T> lr01 = il::lowRank(epsilon, A, s01);
    bie::LowRank<T> lr10 = il::lowRank(epsilon, A, s10);
    bie::LowRank<T> lr11 = il::lowRank(epsilon, A, s11);
    const il::int_t n00 = lr00.A.size(0);
    const il::int_t n10 = lr10.A.size(0);
    const il::int_t n01 = lr00.B.size(0);
    const il::int_t n11 = lr01.B.size(0);
    const il::int_t n0 = n00 + n10;
    const il::int_t n1 = n01 + n11;
    const il::int_t r00 = lr00.A.size(1);
    const il::int_t r01 = lr01.A.size(1);
    const il::int_t r10 = lr10.A.size(1);
    const il::int_t r11 = lr11.A.size(1);

    il::Array2D<T> m00a{n0, r00, 0.0};
    il::Array2D<T> m00b{n1, r00, 0.0};
    il::copy(lr00.A.view(), il::io,
             m00a.Edit(il::Range{0, n00}, il::Range{0, r00}));
    il::copy(lr00.B.view(), il::io,
             m00b.Edit(il::Range{0, n01}, il::Range{0, r00}));
    il::Array2D<T> m01a{n0, r01, 0.0};
    il::Array2D<T> m01b{n1, r01, 0.0};
    il::copy(lr01.A.view(), il::io,
             m01a.Edit(il::Range{0, n00}, il::Range{0, r01}));
    il::copy(lr01.B.view(), il::io,
             m01b.Edit(il::Range{n01, n01 + n11}, il::Range{0, r01}));
    bie::LowRank<T> lr0 = il::lowRankAddition(
        epsilon, 1.0, m00a.view(), m00b.view(), 1.0, m01a.view(), m01b.view());

    il::Array2D<T> m10a{n0, r10, 0.0};
    il::Array2D<T> m10b{n1, r10, 0.0};
    il::copy(lr10.A.view(), il::io,
             m10a.Edit(il::Range{n00, n00 + n10}, il::Range{0, r10}));
    il::copy(lr10.B.view(), il::io,
             m10b.Edit(il::Range{0, n01}, il::Range{0, r10}));
    il::Array2D<T> m11a{n0, r11, 0.0};
    il::Array2D<T> m11b{n1, r11, 0.0};
    il::copy(lr11.A.view(), il::io,
             m11a.Edit(il::Range{n00, n00 + n10}, il::Range{0, r11}));
    il::copy(lr11.B.view(), il::io,
             m11b.Edit(il::Range{n01, n01 + n11}, il::Range{0, r11}));
    bie::LowRank<T> lr1 = il::lowRankAddition(
        epsilon, 1.0, m10a.view(), m10b.view(), 1.0, m11a.view(), m11b.view());

    bie::LowRank<T> ans =
        il::lowRankAddition(epsilon, 1.0, lr0.A.view(), lr0.B.view(), 1.0,
                            lr1.A.view(), lr1.B.view());
    return ans;
  }
}

}  // namespace il