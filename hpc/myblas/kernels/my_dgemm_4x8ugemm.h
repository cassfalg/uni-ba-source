#ifndef HPC_MYBLAS_MY_DGEMM_4X8UGEMM_H
#define HPC_MYBLAS_MY_DGEMM_4X8UGEMM_H 1

#include <hpc/myblas/blocksize.h>
static_assert(hpc::myblas::BlockSize<double>::MR == 4, "MR must be 4 for this micro kernel.");
static_assert(hpc::myblas::BlockSize<double>::NR == 8, "NR must be 8 for this micro kernel.");

extern "C" {
	void my_dgemm_4x8ugemm(
			long k,
			double alpha,
			const double *A,
			const double *B,
			double beta,
			double *C,
			long incRowC,
			long incColC
			);
}

template <typename Index, typename T>
void
ugemm(Index kc, T alpha,
      const T *A, const T *B,
      T beta,
      T *C, Index incRowC, Index incColC)
{
	my_dgemm_4x8ugemm(kc, alpha, A, B, beta, C, incRowC, incColC);
}

#endif // HPC_MYBLAS_MY_DGEMM_4X8UGEMM_H
