#ifndef HPC_REFBLAS_GEMM_H
#define HPC_REFBLAS_GEMM_H 1

#include <type_traits>
#include <hpc/myblas/gescal.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace refblas {

//-----------------------------------------------------------------------------


template <typename Index, typename Alpha,
         typename TA, typename TB,
         typename Beta,
         typename TC>
void
gemm(Index m, Index n, Index k,
     Alpha alpha,
     const TA *A, Index incRowA, Index incColA,
     const TB *B, Index incRowB, Index incColB,
     Beta beta,
     TC *C, Index incRowC, Index incColC)
{

    if (alpha==Alpha(0) || k==0) {
        hpc::myblas::gescal(m, n, beta, C, incRowC, incColC);
        return;
    } else {
        hpc::myblas::gescal(m, n, beta, C, incRowC, incColC);
    }

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            for (Index l=0; l<k; ++l) {
                C[i*incRowC+j*incColC] += alpha*A[i*incRowA+l*incColA]
                                               *B[l*incRowB+j*incColB];
            }
        }
    }
}

template <typename Alpha, typename MA, typename MB, typename Beta, typename MC>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value
                     && hpc::matvec::IsGeMatrix<MB>::value
                     && hpc::matvec::IsGeMatrix<MC>::value,
         void>::type
mm(const Alpha &alpha, const MA &A, const MB &B, const Beta &beta, MC &C)
{
    assert(A.numCols==B.numRows);
    assert(C.numRows>=A.numRows);
    assert(C.numCols>=B.numCols);

    typedef typename std::common_type<typename MA::Index,
                                      typename MB::Index,
                                      typename MC::Index>::type  Index;

    const Index m = C.numRows;
    const Index n = C.numCols;
    const Index k = A.numCols;

    gemm(m, n, k,
                  alpha,
                  A.data, A.incRow, A.incCol,
                  B.data, B.incRow, B.incCol,
                  beta,
                  C.data, C.incRow, C.incCol);
}

} } // namespace refblas, hpc

#endif // HPC_REFBLAS_GEMM_REF_H
