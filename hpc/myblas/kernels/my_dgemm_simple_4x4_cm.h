#ifndef HPC_MYBLAS_MY_DGEMM_SIMPLE_4X4_CM_H
#define HPC_MYBLAS_MY_DGEMM_SIMPLE_4X4_CM_H 1

//#include <hpc/myblas/blocksize.h>
//static_assert(MR == 4, "MR must be 4 for this micro kernel.");
//static_assert(NR == 4, "NR must be 4 for this micro kernel.");

extern "C" {
	void my_dgemm_simple_4x4_cm(
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
	my_dgemm_simple_4x4_cm(kc, alpha, A, B, beta, C, incRowC, incColC);
}

#endif // HPC_MYBLAS_MY_DGEMM_SIMPLE_4X4_CM_H
