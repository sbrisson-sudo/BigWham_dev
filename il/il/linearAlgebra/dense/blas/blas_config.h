//==============================================================================
//
// Copyright 2018 The InsideLoop Authors. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//==============================================================================

#ifndef IL_BLAS_DEFINE_H
#define IL_BLAS_DEFINE_H

#ifdef IL_MKL

#include <mkl_cblas.h>
#define IL_CBLAS_INT MKL_INT
#define IL_CBLAS_LAYOUT CBLAS_LAYOUT
#define IL_CBLAS_PCOMPLEX64 float*
#define IL_CBLAS_PCOMPLEX128 double*
#define IL_CBLAS_PCOMPLEX64_ANS float*
#define IL_CBLAS_PCOMPLEX128_ANS double*

#elif IL_OPENBLAS

#include <cblas.h>
#define IL_CBLAS_INT int
#define IL_CBLAS_LAYOUT CBLAS_ORDER
#define IL_CBLAS_PCOMPLEX64 float*
#define IL_CBLAS_PCOMPLEX128 double*
#define IL_CBLAS_PCOMPLEX64_ANS openblas_complex_float*
#define IL_CBLAS_PCOMPLEX128_ANS openblas_complex_double*

#elif IL_APPLE_ACCELERATE

#include <Accelerate/Accelerate.h>
#define IL_CBLAS_INT int
#define IL_CBLAS_LAYOUT CBLAS_ORDER
#define IL_CBLAS_PCOMPLEX64 float*
#define IL_CBLAS_PCOMPLEX128 double*
#define IL_CBLAS_PCOMPLEX64_ANS float*
#define IL_CBLAS_PCOMPLEX128_ANS double*

// Apple Accelerate uses int for LAPACK operations
typedef int lapack_int;
// Define MKL_INT for compatibility with existing code
typedef int MKL_INT;

// Apple Accelerate AXPBY implementations (not available in Accelerate, so we implement them)
inline void cblas_saxpby(const int N, const float alpha, const float *X, const int incX, const float beta, float *Y, const int incY) {
  cblas_sscal(N, beta, Y, incY);
  cblas_saxpy(N, alpha, X, incX, Y, incY);
}

inline void cblas_daxpby(const int N, const double alpha, const double *X, const int incX, const double beta, double *Y, const int incY) {
  cblas_dscal(N, beta, Y, incY);
  cblas_daxpy(N, alpha, X, incX, Y, incY);
}

inline void cblas_caxpby(const int N, const void *alpha, const void *X, const int incX, const void *beta, void *Y, const int incY) {
  cblas_cscal(N, beta, Y, incY);
  cblas_caxpy(N, alpha, X, incX, Y, incY);
}

inline void cblas_zaxpby(const int N, const void *alpha, const void *X, const int incX, const void *beta, void *Y, const int incY) {
  cblas_zscal(N, beta, Y, incY);
  cblas_zaxpy(N, alpha, X, incX, Y, incY);
}

// Apple Accelerate LAPACK definitions
#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

// Apple Accelerate LAPACKE interface wrappers using existing clapack functions
inline int LAPACKE_dgesv(int layout, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
  __CLPK_integer n_ = n, nrhs_ = nrhs, lda_ = lda, ldb_ = ldb, info = 0;
  dgesv_(&n_, &nrhs_, a, &lda_, ipiv, b, &ldb_, &info);
  return info;
}

inline int LAPACKE_dgbsv(int layout, int n, int kl, int ku, int nrhs, double* ab, int ldab, int* ipiv, double* b, int ldb) {
  __CLPK_integer n_ = n, kl_ = kl, ku_ = ku, nrhs_ = nrhs, ldab_ = ldab, ldb_ = ldb, info = 0;
  dgbsv_(&n_, &kl_, &ku_, &nrhs_, ab, &ldab_, ipiv, b, &ldb_, &info);
  return info;
}

inline int LAPACKE_dgtsv(int layout, int n, int nrhs, double* dl, double* d, double* du, double* b, int ldb) {
  __CLPK_integer n_ = n, nrhs_ = nrhs, ldb_ = ldb, info = 0;
  dgtsv_(&n_, &nrhs_, dl, d, du, b, &ldb_, &info);
  return info;
}

// Additional LAPACKE wrappers needed by the codebase
inline int LAPACKE_dgetrf(int layout, int m, int n, double* a, int lda, int* ipiv) {
  __CLPK_integer m_ = m, n_ = n, lda_ = lda, info = 0;
  dgetrf_(&m_, &n_, a, &lda_, ipiv, &info);
  return info;
}

inline int LAPACKE_dgetrs(int layout, char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb) {
  __CLPK_integer n_ = n, nrhs_ = nrhs, lda_ = lda, ldb_ = ldb, info = 0;
  dgetrs_(&trans, &n_, &nrhs_, const_cast<double*>(a), &lda_, const_cast<int*>(ipiv), b, &ldb_, &info);
  return info;
}

inline int LAPACKE_dgetri(int layout, int n, double* a, int lda, const int* ipiv) {
  __CLPK_integer n_ = n, lda_ = lda, lwork = -1, info = 0;
  double work_query;
  dgetri_(&n_, a, &lda_, const_cast<int*>(ipiv), &work_query, &lwork, &info);
  lwork = static_cast<__CLPK_integer>(work_query);
  double* work = new double[lwork];
  dgetri_(&n_, a, &lda_, const_cast<int*>(ipiv), work, &lwork, &info);
  delete[] work;
  return info;
}

inline int LAPACKE_dgecon(int layout, char norm, int n, const double* a, int lda, double anorm, double* rcond) {
  __CLPK_integer n_ = n, lda_ = lda, info = 0;
  double* work = new double[4*n];
  __CLPK_integer* iwork = new __CLPK_integer[n];
  dgecon_(&norm, &n_, const_cast<double*>(a), &lda_, &anorm, rcond, work, iwork, &info);
  delete[] work;
  delete[] iwork;
  return info;
}

// Complex number type definitions for Apple Accelerate
typedef __CLPK_doublecomplex lapack_complex_double;
typedef __CLPK_complex lapack_complex_float;

// Additional critical LAPACK functions
inline int LAPACKE_dpotrf(int layout, char uplo, int n, double* a, int lda) {
  __CLPK_integer n_ = n, lda_ = lda, info = 0;
  dpotrf_(&uplo, &n_, a, &lda_, &info);
  return info;
}

inline int LAPACKE_dpotrs(int layout, char uplo, int n, int nrhs, const double* a, int lda, double* b, int ldb) {
  __CLPK_integer n_ = n, nrhs_ = nrhs, lda_ = lda, ldb_ = ldb, info = 0;
  dpotrs_(&uplo, &n_, &nrhs_, const_cast<double*>(a), &lda_, b, &ldb_, &info);
  return info;
}

inline int LAPACKE_dpotri(int layout, char uplo, int n, double* a, int lda) {
  __CLPK_integer n_ = n, lda_ = lda, info = 0;
  dpotri_(&uplo, &n_, a, &lda_, &info);
  return info;
}

inline int LAPACKE_dpocon(int layout, char uplo, int n, const double* a, int lda, double anorm, double* rcond) {
  __CLPK_integer n_ = n, lda_ = lda, info = 0;
  double* work = new double[3*n];
  __CLPK_integer* iwork = new __CLPK_integer[n];
  dpocon_(&uplo, &n_, const_cast<double*>(a), &lda_, &anorm, rcond, work, iwork, &info);
  delete[] work;
  delete[] iwork;
  return info;
}

// Proper implementations for Apple Accelerate
inline int LAPACKE_dgesvd(int layout, char jobu, char jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double* superb) {
  __CLPK_integer m_ = m, n_ = n, lda_ = lda, ldu_ = ldu, ldvt_ = ldvt, info = 0;
  __CLPK_integer lwork = -1;
  double work_query;
  
  // Query optimal workspace size
  dgesvd_(&jobu, &jobvt, &m_, &n_, a, &lda_, s, u, &ldu_, vt, &ldvt_, &work_query, &lwork, &info);
  
  if (info != 0) return info;
  
  lwork = static_cast<__CLPK_integer>(work_query);
  double* work = new double[lwork];
  
  // Perform SVD
  dgesvd_(&jobu, &jobvt, &m_, &n_, a, &lda_, s, u, &ldu_, vt, &ldvt_, work, &lwork, &info);
  
  delete[] work;
  return info;
}

inline int LAPACKE_dgebrd(int layout, int m, int n, double* a, int lda, double* d, double* e, double* tauq, double* taup) {
  __CLPK_integer m_ = m, n_ = n, lda_ = lda, info = 0;
  __CLPK_integer lwork = -1;
  double work_query;
  
  // Query optimal workspace size
  dgebrd_(&m_, &n_, a, &lda_, d, e, tauq, taup, &work_query, &lwork, &info);
  
  if (info != 0) return info;
  
  lwork = static_cast<__CLPK_integer>(work_query);
  double* work = new double[lwork];
  
  // Perform bidiagonal reduction
  dgebrd_(&m_, &n_, a, &lda_, d, e, tauq, taup, work, &lwork, &info);
  
  delete[] work;
  return info;
}

inline int LAPACKE_dbdsqr(int layout, char uplo, int n, int ncvt, int nru, int ncc, double* d, double* e, double* vt, int ldvt, double* u, int ldu, double* c, int ldc) {
  __CLPK_integer n_ = n, ncvt_ = ncvt, nru_ = nru, ncc_ = ncc;
  __CLPK_integer ldvt_ = ldvt, ldu_ = ldu, ldc_ = ldc, info = 0;
  double* work = new double[4*n];
  
  dbdsqr_(&uplo, &n_, &ncvt_, &nru_, &ncc_, d, e, vt, &ldvt_, u, &ldu_, c, &ldc_, work, &info);
  
  delete[] work;
  return info;
}

inline int LAPACKE_dpptrf(int layout, char uplo, int n, double* ap) {
  __CLPK_integer n_ = n, info = 0;
  dpptrf_(&uplo, &n_, ap, &info);
  return info;
}

inline int LAPACKE_dgehrd(int layout, int n, int ilo, int ihi, double* a, int lda, double* tau) {
  __CLPK_integer n_ = n, ilo_ = ilo, ihi_ = ihi, lda_ = lda, info = 0;
  __CLPK_integer lwork = -1;
  double work_query;
  
  // Query optimal workspace size
  dgehrd_(&n_, &ilo_, &ihi_, a, &lda_, tau, &work_query, &lwork, &info);
  
  if (info != 0) return info;
  
  lwork = static_cast<__CLPK_integer>(work_query);
  double* work = new double[lwork];
  
  // Perform Hessenberg reduction
  dgehrd_(&n_, &ilo_, &ihi_, a, &lda_, tau, work, &lwork, &info);
  
  delete[] work;
  return info;
}

inline int LAPACKE_dhseqr(int layout, char job, char compz, int n, int ilo, int ihi, double* h, int ldh, double* wr, double* wi, double* z, int ldz) {
  __CLPK_integer n_ = n, ilo_ = ilo, ihi_ = ihi, ldh_ = ldh, ldz_ = ldz, info = 0;
  __CLPK_integer lwork = -1;
  double work_query;
  
  // Query optimal workspace size
  dhseqr_(&job, &compz, &n_, &ilo_, &ihi_, h, &ldh_, wr, wi, z, &ldz_, &work_query, &lwork, &info);
  
  if (info != 0) return info;
  
  lwork = static_cast<__CLPK_integer>(work_query);
  double* work = new double[lwork];
  
  // Perform eigenvalue computation
  dhseqr_(&job, &compz, &n_, &ilo_, &ihi_, h, &ldh_, wr, wi, z, &ldz_, work, &lwork, &info);
  
  delete[] work;
  return info;
}

// Complex number implementations
inline int LAPACKE_zgetrf(int layout, int m, int n, lapack_complex_double* a, int lda, int* ipiv) {
  __CLPK_integer m_ = m, n_ = n, lda_ = lda, info = 0;
  zgetrf_(&m_, &n_, reinterpret_cast<__CLPK_doublecomplex*>(a), &lda_, ipiv, &info);
  return info;
}

inline int LAPACKE_zgetri(int layout, int n, lapack_complex_double* a, int lda, const int* ipiv) {
  __CLPK_integer n_ = n, lda_ = lda, lwork = -1, info = 0;
  __CLPK_doublecomplex work_query;
  
  // Query optimal workspace size
  zgetri_(&n_, reinterpret_cast<__CLPK_doublecomplex*>(a), &lda_, const_cast<int*>(ipiv), &work_query, &lwork, &info);
  
  if (info != 0) return info;
  
  lwork = static_cast<__CLPK_integer>(work_query.r);
  __CLPK_doublecomplex* work = new __CLPK_doublecomplex[lwork];
  
  // Perform matrix inverse
  zgetri_(&n_, reinterpret_cast<__CLPK_doublecomplex*>(a), &lda_, const_cast<int*>(ipiv), work, &lwork, &info);
  
  delete[] work;
  return info;
}

#endif

#endif  // IL_BLAS_DEFINE_H
