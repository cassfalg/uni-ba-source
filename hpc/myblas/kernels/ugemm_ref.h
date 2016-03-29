#ifndef HPC_MYBLAS_KERNELS_UGEMM_REF_H
#define HPC_MYBLAS_KERNELS_UGEMM_REF_H 1

#include <hpc/myblas/blocksize.h>

namespace hpc { namespace myblas {

// ugemm with beta == 0
  template <typename Index, typename T>
  void
  ugemm_0(Index kc, T alpha,
        const T *A, const T *B,
        T *C, Index incRowC, Index incColC)
  {
      const Index MR = BlockSize<T>::MR;
      const Index NR = BlockSize<T>::NR;
      T P[MR*NR];

      for (Index l=0; l<MR*NR; ++l) {
          P[l] = 0;
      }
      for (Index l=0; l<kc; ++l) {
          for (Index j=0; j<NR; ++j) {
              for (Index i=0; i<MR; ++i) {
                  P[i+j*MR] += A[i+l*MR]*B[l*NR+j];
              }
          }
      }
      for (Index j=0; j<NR; ++j) {
          for (Index i=0; i<MR; ++i) {
              C[i*incRowC+j*incColC] += alpha*P[i+j*MR];
          }
      }
  }

template <typename Index, typename T>
void
ugemm(Index kc, T alpha,
      const T *A, const T *B,
      T beta,
      T *C, Index incRowC, Index incColC)
{
    const Index MR = BlockSize<T>::MR;
    const Index NR = BlockSize<T>::NR;
    T P[MR*NR];

    for (Index l=0; l<MR*NR; ++l) {
        P[l] = 0;
    }
    for (Index l=0; l<kc; ++l) {
        for (Index j=0; j<NR; ++j) {
            for (Index i=0; i<MR; ++i) {
                P[i+j*MR] += A[i+l*MR]*B[l*NR+j];
            }
        }
    }
    for (Index j=0; j<NR; ++j) {
        for (Index i=0; i<MR; ++i) {
            C[i*incRowC+j*incColC] *= beta;
            C[i*incRowC+j*incColC] += alpha*P[i+j*MR];
        }
    }
}

} } // namespace mvblas, hpc

#endif // HPC_MYBLAS_KERNELS_UGEMM_REF_H
