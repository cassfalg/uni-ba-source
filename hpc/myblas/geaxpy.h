#ifndef HPC_MYBLAS_GEAXPY_H
#define HPC_MYBLAS_GEAXPY_H 1

namespace hpc { namespace myblas {

template <typename Index, typename Alpha, typename TX, typename TY>
void
geaxpy(Index m, Index n,
       const Alpha &alpha,
       const TX *X, Index incRowX, Index incColX,
       TY       *Y, Index incRowY, Index incColY)
{
    if (incRowY<incColY) {
        for (Index j=0; j<n; ++j) {
            for (Index i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    } else {
        for (Index i=0; i<m; ++i) {
            for (Index j=0; j<n; ++j) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    }
}

} } // namespace myblas, hpc

#endif // HPC_MYBLAS_GEAXPY_H
