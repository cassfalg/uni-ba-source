#include <chrono>
#include <complex>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/myblas/gemm.h>
#include <hpc/refblas/gemm.h>
#include <hpc/myblas/blocksize.h>
#include <hpc/matvec/random_init.h>
#include <hpc/util/walltime.h>

//    1: single precision
//    2: double precision
#ifndef PRECISION
#define PRECISION 2
#endif


int
main()
{
	using namespace hpc::matvec;
#if PRECISION == 1
    typedef float ElementType;
	#ifndef SIZE
	#define SIZE 1400u
	#endif
#elif PRECISION == 2
    typedef double ElementType;
	#ifndef SIZE
	#define SIZE 1000u
	#endif
#endif

    typedef size_t Index;
    typedef GeMatrix<ElementType, Index> MT;

    MT A(SIZE, SIZE, RowMajor);
    MT B(SIZE, SIZE, ColMajor);
    MT C(SIZE, SIZE, ColMajor);

    randomInit(A);
    randomInit(B);
    randomInit(C);
    ElementType alpha(1.5);
    ElementType beta(2.5);

    hpc::util::WallTime<double> wallTime;

#if PRECISION == 1
    printf("SP: %5d %5d %5d %5d %5d ", S_BLOCKSIZE_MC, S_BLOCKSIZE_KC,
	   S_BLOCKSIZE_NC, S_BLOCKSIZE_MR, S_BLOCKSIZE_NR);
#elif PRECISION == 2
    printf("DP: %5d %5d %5d %5d %5d ", D_BLOCKSIZE_MC, D_BLOCKSIZE_KC,
       D_BLOCKSIZE_NC, D_BLOCKSIZE_MR, D_BLOCKSIZE_NR);
#endif

    wallTime.tic();
    hpc::myblas::gemm(Index(SIZE), Index(SIZE), Index(SIZE), alpha,
                  A.data, A.incRow, A.incCol,
                  B.data, B.incRow, B.incCol,
                  beta,
                  C.data, C.incRow, C.incCol);
    double t = wallTime.toc();

    printf("    %20.4lf %9.2lf ", t, 2.*SIZE/1000*SIZE/1000*SIZE/t);
    printf("\n");

}
