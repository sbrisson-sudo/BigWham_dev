//
// Created by peruzzo on 03.05.21.
//

#ifndef IL_SVD_H
#define IL_SVD_H

#ifdef IL_MKL
#include <mkl_lapacke.h>
#elif IL_OPENBLAS
#include <OpenBLAS/lapacke.h>
#endif

namespace il {
inline void svdDecomposition(il::io_t, il::Array2DEdit<double> A, il::Array2DEdit<double> A, ) {

    const int layout = LAPACK_COL_MAJOR;
    //const int layout = LAPACK_ROW_MAJOR;

    const lapack_int m = static_cast<lapack_int>(A.size(0));
    const lapack_int n = static_cast<lapack_int>(A.size(1));

    const char jobu = 'A';
    const char jobvt = 'A';

    //const lapack_int lda = static_cast<lapack_int>(A.stride(1)); // to judge if A.stride(1) is min(m,n)
    const il::int_t min_mn = m < n ? m : n;
    const lapack_int lda = static_cast<lapack_int>(min_mn);

    double superb[min_mn-1];
    /* Local arrays */
    double s[N], u[LDU*M], vt[LDVT*N];

    /* Compute SVD */
    const lapack_int lapack_error = = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, jobu, jobvt, m, n, A.Data(), lda, s, u, ldu, vt, ldvt, superb );

    IL_EXPECT_FAST(lapack_error == 0);
}

}  // namespace il
#endif //IL_SVD_H
