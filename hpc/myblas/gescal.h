#ifndef HPC_MYBLAS_GESCAL_H
#define HPC_MYBLAS_GESCAL_H 1

namespace hpc { namespace myblas {

template <typename Index, typename Alpha, typename TX>
void
gescal(Index m, Index n,
       const Alpha &alpha,
       TX *X, Index incRowX, Index incColX)
{
    if (alpha!=Alpha(1)) {
        if (incRowX<incColX) {
            for (Index j=0; j<n; ++j) {
                for (Index i=0; i<m; ++i) {
                    X[i*incRowX+j*incColX] *= alpha;
                }
            }
        } else {
            for (Index i=0; i<m; ++i) {
                for (Index j=0; j<n; ++j) {
                    X[i*incRowX+j*incColX] *= alpha;
                }
            }
        }
    }
}

} } // namespace myblas, hpc

#endif // HPC_MYBLAS_GESCAL_H
