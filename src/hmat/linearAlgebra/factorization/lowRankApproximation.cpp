#include <il/Array.h>

#include <il/blas.h>
#include <il/linearAlgebra/Matrix.h>
#include <il/linearAlgebra/dense/factorization/qrDecomposition.h>

#include <hmat/linearAlgebra/factorization/lowRankApproximation.h>

#ifdef IL_MKL
#include <mkl_lapacke.h>
#elif IL_OPENBLAS
#include <OpenBLAS/lapacke.h>
#endif

namespace il {

bie::LowRank<double> lowRank(double epsilon, il::Array2DView<double> A) {
  const il::int_t n0 = A.size(0);
  const il::int_t n1 = A.size(1);
  il::Array2D<double> M{n0, n1};
  il::copy(A, il::io, M.Edit());

  il::Array2D<double> U{n0, n0};
  il::Array2D<double> V{n1, n1};
  il::Array<double> svalue{il::min(n0, n1)};
  il::Array<double> superb{il::min(n0, n1) - 1};
  {
    const int layout = LAPACK_COL_MAJOR;
    const char jobu = 'A';
    const char jobvt = 'A';
    const lapack_int m = static_cast<lapack_int>(n0);
    const lapack_int n = static_cast<lapack_int>(n1);
    const lapack_int lda = static_cast<lapack_int>(M.stride(1));
    const lapack_int ldu = static_cast<lapack_int>(U.stride(1));
    const lapack_int ldvt = static_cast<lapack_int>(V.stride(1));
    const lapack_int lapack_error =
        LAPACKE_dgesvd(layout, jobu, jobvt, m, n, M.Data(), lda, svalue.Data(),
                       U.Data(), ldu, V.Data(), ldvt, superb.Data());
    IL_EXPECT_FAST(lapack_error == 0);
  }

  il::Array2D<double> a{n0, 0, 0.0};
  il::Array2D<double> b{n1, 0, 0.0};
  {
    il::int_t k = 0;
    while (k < n0 && k < n1 && svalue[k] >= epsilon * svalue[0]) {
      a.Resize(a.size(0), a.size(1) + 1, 0.0);
      b.Resize(b.size(0), b.size(1) + 1, 0.0);
      for (il::int_t i = 0; i < n0; ++i) {
        a(i, k) = std::sqrt(svalue[k]) * U(i, k);
      }
      for (il::int_t i = 0; i < n1; ++i) {
        b(i, k) = std::sqrt(svalue[k]) * V(k, i);
      }
      ++k;
    }
  }

  IL_EXPECT_FAST(a.size(1) > 0);

  bie::LowRank<double> ans{};
  ans.A = std::move(a);
  ans.B = std::move(b);

  return ans;
}

bie::LowRank<std::complex<double>> lowRank(
    double epsilon, il::Array2DView<std::complex<double>> A) {
  const il::int_t n0 = A.size(0);
  const il::int_t n1 = A.size(1);
  il::Array2D<std::complex<double>> M{n0, n1};
  il::copy(A, il::io, M.Edit());

  il::Array2D<std::complex<double>> U{n0, n0};
  il::Array2D<std::complex<double>> V{n1, n1};
  il::Array<double> svalue{il::min(n0, n1)};
  il::Array<double> superb{il::min(n0, n1) - 1};
  {
    const int layout = LAPACK_COL_MAJOR;
    const char jobu = 'A';
    const char jobvt = 'A';
    const lapack_int m = static_cast<lapack_int>(n0);
    const lapack_int n = static_cast<lapack_int>(n1);
    const lapack_int lda = static_cast<lapack_int>(M.stride(1));
    const lapack_int ldu = static_cast<lapack_int>(U.stride(1));
    const lapack_int ldvt = static_cast<lapack_int>(V.stride(1));
    const lapack_int lapack_error = LAPACKE_zgesvd(
        layout, jobu, jobvt, m, n,
        reinterpret_cast<lapack_complex_double*>(M.Data()), lda, svalue.Data(),
        reinterpret_cast<lapack_complex_double*>(U.Data()), ldu,
        reinterpret_cast<lapack_complex_double*>(V.Data()), ldvt,
        superb.Data());
    IL_EXPECT_FAST(lapack_error == 0);
  }

  il::Array2D<std::complex<double>> a{n0, 0, 0.0};
  il::Array2D<std::complex<double>> b{n1, 0, 0.0};
  {
    il::int_t k = 0;
    while (k < n0 && k < n1 && svalue[k] >= epsilon * svalue[0]) {
      a.Resize(a.size(0), a.size(1) + 1, 0.0);
      b.Resize(b.size(0), b.size(1) + 1, 0.0);
      for (il::int_t i = 0; i < n0; ++i) {
        a(i, k) = std::sqrt(svalue[k]) * U(i, k);
      }
      for (il::int_t i = 0; i < n1; ++i) {
        b(i, k) = std::sqrt(svalue[k]) * V(k, i);
      }
      ++k;
    }
  }

  IL_EXPECT_FAST(a.size(1) > 0);

  bie::LowRank<std::complex<double>> ans{};
  ans.A = std::move(a);
  ans.B = std::move(b);

  return ans;
}

bie::LowRank<double> lowRankAddition(double epsilon, double alpha,
                                    il::Array2DView<double> aa,
                                    il::Array2DView<double> ab, double beta,
                                    il::Array2DView<double> ba,
                                    il::Array2DView<double> bb) {
  IL_EXPECT_MEDIUM(aa.size(1) == ab.size(1));
  IL_EXPECT_MEDIUM(ba.size(1) == bb.size(1));
  IL_EXPECT_MEDIUM(aa.size(0) == ba.size(0));
  IL_EXPECT_MEDIUM(ab.size(0) == bb.size(0));

  const il::int_t na = aa.size(0);
  const il::int_t nb = ab.size(0);
  const il::int_t ra = aa.size(1);
  const il::int_t rb = ba.size(1);

  il::Array2D<double> la{na, ra + rb};
  il::blas(alpha, aa, 0.0, il::io, la.Edit(il::Range{0, na}, il::Range{0, ra}));
  il::blas(beta, ba, 0.0, il::io,
           la.Edit(il::Range{0, na}, il::Range{ra, ra + rb}));
  const il::int_t pa = il::min(na, ra + rb);
  il::Array<double> taua{pa};
  il::qrDecomposition(il::io, taua.Edit(), la.Edit());

  il::Array2D<double> lb{nb, ra + rb};
  il::copy(ab, il::io, lb.Edit(il::Range{0, nb}, il::Range{0, ra}));
  il::copy(bb, il::io, lb.Edit(il::Range{0, nb}, il::Range{ra, ra + rb}));
  const il::int_t pb = il::min(nb, ra + rb);
  il::Array<double> taub{pb};
  il::qrDecomposition(il::io, taub.Edit(), lb.Edit());

  il::Array2D<double> P{pa, pb};
  for (il::int_t i1 = 0; i1 < pb; ++i1) {
    for (il::int_t i0 = 0; i0 < pa; ++i0) {
      double sum = 0.0;
      for (il::int_t k = il::max(i0, i1); k < ra + rb; ++k) {
        sum += la(i0, k) * lb(i1, k);
      }
      P(i0, i1) = sum;
    }
  }

  // Begin SVD compression
  il::Array2D<double> U{pa, pa};
  il::Array2D<double> V{pb, pb};
  il::Array<double> svalue{il::min(pa, pb)};
  il::Array<double> superb{il::min(pa, pb) - 1};
  {
    const int layout = LAPACK_COL_MAJOR;
    const char jobu = 'A';
    const char jobvt = 'A';
    const lapack_int m = static_cast<lapack_int>(P.size(0));
    const lapack_int n = static_cast<lapack_int>(P.size(1));
    const lapack_int lda = static_cast<lapack_int>(P.stride(1));
    const lapack_int ldu = static_cast<lapack_int>(U.stride(1));
    const lapack_int ldvt = static_cast<lapack_int>(V.stride(1));
    const lapack_int lapack_error =
        LAPACKE_dgesvd(layout, jobu, jobvt, m, n, P.Data(), lda, svalue.Data(),
                       U.Data(), ldu, V.Data(), ldvt, superb.Data());
    IL_EXPECT_FAST(lapack_error == 0);
  }

  il::Array2D<double> a{na, 0, 0.0};
  il::Array2D<double> b{nb, 0, 0.0};
  {
    il::int_t k = 0;
    while (k < pa && k < pb && svalue[k] >= epsilon * svalue[0]) {
      a.Resize(a.size(0), a.size(1) + 1, 0.0);
      b.Resize(b.size(0), b.size(1) + 1, 0.0);
      for (il::int_t i = 0; i < pa; ++i) {
        a(i, k) = std::sqrt(svalue[k]) * U(i, k);
      }
      for (il::int_t i = 0; i < pb; ++i) {
        b(i, k) = std::sqrt(svalue[k]) * V(k, i);
      }
      ++k;
    }
  }
  // End SVD compression

  // a <- Qa.a
  {
    const int layout = LAPACK_COL_MAJOR;
    const char side = 'L';
    const char trans = 'N';
    const lapack_int m = static_cast<lapack_int>(a.size(0));
    const lapack_int n = static_cast<lapack_int>(a.size(1));
    const lapack_int k = static_cast<lapack_int>(pa);
    const lapack_int ldla = static_cast<lapack_int>(la.stride(1));
    const lapack_int lda = static_cast<lapack_int>(a.stride(1));
    const lapack_int lapack_error =
        LAPACKE_dormqr(layout, side, trans, m, n, k, la.data(), ldla,
                       taua.data(), a.Data(), lda);
    IL_EXPECT_FAST(lapack_error == 0);
  }
  // b <- Qb.b
  {
    const int layout = LAPACK_COL_MAJOR;
    const char side = 'L';
    const char trans = 'N';
    const lapack_int m = static_cast<lapack_int>(b.size(0));
    const lapack_int n = static_cast<lapack_int>(b.size(1));
    const lapack_int k = static_cast<lapack_int>(pb);
    const lapack_int ldlb = static_cast<lapack_int>(lb.stride(1));
    const lapack_int ldb = static_cast<lapack_int>(b.stride(1));
    const lapack_int lapack_error =
        LAPACKE_dormqr(layout, side, trans, m, n, k, lb.data(), ldlb,
                       taub.data(), b.Data(), ldb);
    IL_EXPECT_FAST(lapack_error == 0);
  }

  bie::LowRank<double> ans{};
  ans.A = std::move(a);
  ans.B = std::move(b);

  return ans;
}

bie::LowRank<std::complex<double>> lowRankAddition(
    double epsilon, std::complex<double> alpha,
    il::Array2DView<std::complex<double>> aa,
    il::Array2DView<std::complex<double>> ab, std::complex<double> beta,
    il::Array2DView<std::complex<double>> ba,
    il::Array2DView<std::complex<double>> bb) {
  IL_EXPECT_MEDIUM(aa.size(1) == ab.size(1));
  IL_EXPECT_MEDIUM(ba.size(1) == bb.size(1));
  IL_EXPECT_MEDIUM(aa.size(0) == ba.size(0));
  IL_EXPECT_MEDIUM(ab.size(0) == bb.size(0));

  const il::int_t na = aa.size(0);
  const il::int_t nb = ab.size(0);
  const il::int_t ra = aa.size(1);
  const il::int_t rb = ba.size(1);

  il::Array2D<std::complex<double>> la{na, ra + rb};
  il::blas(alpha, aa, std::complex<double>{0.0}, il::io,
           la.Edit(il::Range{0, na}, il::Range{0, ra}));
  il::blas(beta, ba, std::complex<double>{0.0}, il::io,
           la.Edit(il::Range{0, na}, il::Range{ra, ra + rb}));
  const il::int_t pa = il::min(na, ra + rb);
  il::Array<std::complex<double>> taua{pa};
  il::qrDecomposition(il::io, taua.Edit(), la.Edit());

  il::Array2D<std::complex<double>> lb{nb, ra + rb};
  il::copy(ab, il::io, lb.Edit(il::Range{0, nb}, il::Range{0, ra}));
  il::copy(bb, il::io, lb.Edit(il::Range{0, nb}, il::Range{ra, ra + rb}));
  const il::int_t pb = il::min(nb, ra + rb);
  il::Array<std::complex<double>> taub{pb};
  il::qrDecomposition(il::io, taub.Edit(), lb.Edit());

  il::Array2D<std::complex<double>> P{pa, pb};
  for (il::int_t i1 = 0; i1 < pb; ++i1) {
    for (il::int_t i0 = 0; i0 < pa; ++i0) {
      std::complex<double> sum = 0.0;
      for (il::int_t k = il::max(i0, i1); k < ra + rb; ++k) {
        sum += la(i0, k) * lb(i1, k);
      }
      P(i0, i1) = sum;
    }
  }

  // Begin SVD compression
  il::Array2D<std::complex<double>> U{pa, pa};
  il::Array2D<std::complex<double>> V{pb, pb};
  il::Array<double> svalue{il::min(pa, pb)};
  il::Array<double> superb{il::min(pa, pb) - 1};
  {
    const int layout = LAPACK_COL_MAJOR;
    const char jobu = 'A';
    const char jobvt = 'A';
    const lapack_int m = static_cast<lapack_int>(P.size(0));
    const lapack_int n = static_cast<lapack_int>(P.size(1));
    const lapack_int lda = static_cast<lapack_int>(P.stride(1));
    const lapack_int ldu = static_cast<lapack_int>(U.stride(1));
    const lapack_int ldvt = static_cast<lapack_int>(V.stride(1));
    const lapack_int lapack_error = LAPACKE_zgesvd(
        layout, jobu, jobvt, m, n,
        reinterpret_cast<lapack_complex_double*>(P.Data()), lda, svalue.Data(),
        reinterpret_cast<lapack_complex_double*>(U.Data()), ldu,
        reinterpret_cast<lapack_complex_double*>(V.Data()), ldvt,
        superb.Data());
    IL_EXPECT_FAST(lapack_error == 0);
  }

  il::Array2D<std::complex<double>> a{na, 0, 0.0};
  il::Array2D<std::complex<double>> b{nb, 0, 0.0};
  {
    il::int_t k = 0;
    while (k < pa && k < pb && svalue[k] >= epsilon * svalue[0]) {
      a.Resize(a.size(0), a.size(1) + 1, 0.0);
      b.Resize(b.size(0), b.size(1) + 1, 0.0);
      for (il::int_t i = 0; i < pa; ++i) {
        a(i, k) = std::sqrt(svalue[k]) * U(i, k);
      }
      for (il::int_t i = 0; i < pb; ++i) {
        b(i, k) = std::sqrt(svalue[k]) * V(k, i);
      }
      ++k;
    }
  }
  // End SVD compression

  // a <- Qa.a
  {
    const int layout = LAPACK_COL_MAJOR;
    const char side = 'L';
    const char trans = 'N';
    const lapack_int m = static_cast<lapack_int>(a.size(0));
    const lapack_int n = static_cast<lapack_int>(a.size(1));
    const lapack_int k = static_cast<lapack_int>(pa);
    const lapack_int ldla = static_cast<lapack_int>(la.stride(1));
    const lapack_int lda = static_cast<lapack_int>(a.stride(1));
    const lapack_int lapack_error = LAPACKE_zunmqr(
        layout, side, trans, m, n, k,
        reinterpret_cast<const lapack_complex_double*>(la.data()), ldla,
        reinterpret_cast<const lapack_complex_double*>(taua.data()),
        reinterpret_cast<lapack_complex_double*>(a.Data()), lda);
    IL_EXPECT_FAST(lapack_error == 0);
  }
  // b <- Qb.b
  {
    const int layout = LAPACK_COL_MAJOR;
    const char side = 'L';
    const char trans = 'N';
    const lapack_int m = static_cast<lapack_int>(b.size(0));
    const lapack_int n = static_cast<lapack_int>(b.size(1));
    const lapack_int k = static_cast<lapack_int>(pb);
    const lapack_int ldlb = static_cast<lapack_int>(lb.stride(1));
    const lapack_int ldb = static_cast<lapack_int>(b.stride(1));
    const lapack_int lapack_error = LAPACKE_zunmqr(
        layout, side, trans, m, n, k,
        reinterpret_cast<const lapack_complex_double*>(lb.data()), ldlb,
        reinterpret_cast<const lapack_complex_double*>(taub.data()),
        reinterpret_cast<lapack_complex_double*>(b.Data()), ldb);
    IL_EXPECT_FAST(lapack_error == 0);
  }

  bie::LowRank<std::complex<double>> ans{};
  ans.A = std::move(a);
  ans.B = std::move(b);

  return ans;
};

}  // namespace il
