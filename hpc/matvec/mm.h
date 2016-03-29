#ifndef HPC_MATVEC_MM_H
#define HPC_MATVEC_MM_H 1

#include <cassert>
#include <type_traits>
#include <hpc/myblas/gemm.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename Alpha, typename MA, typename MB, typename Beta, typename MC>
typename std::enable_if<IsGeMatrix<MA>::value
                     && IsGeMatrix<MB>::value
                     && IsGeMatrix<MC>::value,
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

    myblas::gemm(m, n, k,
                  alpha,
                  A.data, A.incRow, A.incCol,
                  B.data, B.incRow, B.incCol,
                  beta,
                  C.data, C.incRow, C.incCol);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_MM_H
