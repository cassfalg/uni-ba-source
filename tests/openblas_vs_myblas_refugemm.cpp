#include <complex>
#include <cstdio>
#include <hpc/myblas/gecopy.h>
#include <hpc/myblas/gemm.h>
#include <hpc/matvec/residuum_gemm.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/random_init.h>
#include <hpc/util/walltime.h>
#include "openblas/cblas.h"

//    1: single precision
//    2: double precision
#ifndef PRECISION
#define PRECISION 2
#endif

#ifndef MIN_SIZE
#define MIN_SIZE 200u
#endif
#ifndef GRANULARITY
#define GRANULARITY 200u
#endif

int
main() {
	openblas_set_num_threads(1);
    using namespace hpc::matvec;

#if PRECISION == 1
    typedef float ElementType;
	#ifndef MAX_SIZE
	#define MAX_SIZE 10000u // uses about 1,5 GiB of memory
	#endif
#elif PRECISION == 2
    typedef double ElementType;
	#ifndef MAX_SIZE
	#define MAX_SIZE 7000u // uses about 1,5 GiB of memory
	#endif
#endif

    typedef CBLAS_INDEX Index; // openblas uses int? size_t?
    typedef GeMatrix<ElementType, Index> MT;

    MT A(MAX_SIZE, MAX_SIZE);
    MT B(MAX_SIZE, MAX_SIZE);
    MT C1(MAX_SIZE, MAX_SIZE);
    MT C2(MAX_SIZE, MAX_SIZE);

    randomInit(A);
    randomInit(B);
    randomInit(C1);
    ElementType alpha(1.5);
    ElementType beta(2.5);

    hpc::myblas::gecopy<Index, ElementType, ElementType>(MAX_SIZE, MAX_SIZE,
                         C1.data, C1.incRow, C1.incCol,
                         C2.data, C2.incRow, C2.incCol);

    // Header for benchmark
    printf("%5s ", "m=n=k");
#if PRECISION == 1
    printf("%20s %9s", "openblas_sgemm: t", "MFLOPS");
    printf("%20s %9s %9s", "blocked SGEMM: t", "MFLOPS", "diff\n");
#elif PRECISION == 2
    printf("%20s %9s", "openblas_dgemm: t", "MFLOPS");
    printf("%20s %9s %9s", "blocked DGEMM: t", "MFLOPS", "diff\n");
#endif


    hpc::util::WallTime<double> wallTime;

    for (Index i=MIN_SIZE; i<=MAX_SIZE; i+=GRANULARITY) {
        printf("%5lu ", i);

        wallTime.tic();
#if PRECISION == 1
    cblas_sgemm(
#elif PRECISION == 2
    cblas_dgemm(
#endif
		CblasColMajor, CblasNoTrans, CblasNoTrans,
	    i, i, i, alpha, A.data, A.incCol,
	    B.data, B.incCol,
	    beta, C1.data, C1.incCol);
        double t = wallTime.toc();
        printf("%20.4lf %9.2lf", t, 2.*i/1000*i/1000*i/t);

        wallTime.tic();
        hpc::myblas::gemm(i, i, i, alpha,
                      A.data, A.incRow, A.incCol,
                      B.data, B.incRow, B.incCol,
                      beta,
                      C2.data, C2.incRow, C2.incCol);
        t = wallTime.toc();
        double res = hpc::matvec::residuum_gemm(i, i, i, alpha,
                              A.data, A.incRow, A.incCol,
                              B.data, B.incRow, B.incCol,
                              beta,
                              C1.data, C1.incRow, C1.incCol,
                              C2.data, C2.incRow, C2.incCol);

        printf("%20.4lf %9.2lf %9.1e", t, 2.*i/1000*i/1000*i/t, res);
        printf("\n");
    }
}
