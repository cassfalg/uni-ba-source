#include <complex>
#include <cstdio>
#include <hpc/myblas/gecopy.h>
#include <hpc/myblas/gemm.h>
#include <hpc/matvec/residuum_gemm.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/random_init.h>
#include <hpc/util/walltime.h>
extern "C" {
  #include "atlas/cblas.h"
}

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

#if PRECISION == 1
    typedef float ElementType;
	#ifndef MAX_SIZE
	#define MAX_SIZE 9000u // uses about 1,5 GiB of memory, about all we have
	#endif
#elif PRECISION == 2
    typedef double ElementType;
	#ifndef MAX_SIZE
	#define MAX_SIZE 6000u // uses about 1,5 GiB of memory, about all we have
	#endif
#endif

using namespace hpc::matvec;
typedef size_t Index;
typedef GeMatrix<ElementType, Index> MT;

Index granularity(Index i) {
	if (GRANULARITY != 200) {
		return GRANULARITY; // user defined, do not touch
	} else { // use primes, bigger ones with increasing i
		if (i < 260) {
			return 1;
		} else if (i < 800) {
			return 7;
		} else if (i < 2000) {
			return 41;
		} else {
			return 211;
		}
	}
}

int
main() {

    MT A(MAX_SIZE, MAX_SIZE, ColMajor);
    MT B(MAX_SIZE, MAX_SIZE, ColMajor);
    MT C(MAX_SIZE, MAX_SIZE, ColMajor);
    MT C1(MAX_SIZE, MAX_SIZE, ColMajor);
    MT C2(MAX_SIZE, MAX_SIZE, ColMajor);

    randomInit(A);
    randomInit(B);
    randomInit(C);
    ElementType alpha(1.5);
    ElementType beta(2.5);

    // Header for benchmark
    printf("%5s ", "m=n=k");
#if PRECISION == 1
    printf("%20s %9s", "atlas_sgemm: t", "MFLOPS");
    printf("%20s %9s %9s %9s %9s\n", "blocked SGEMM: t", "MFLOPS"
    		, "residuum", "res_alt", "sup-diff");
#elif PRECISION == 2
    printf("%20s %9s", "atlas_dgemm: t", "MFLOPS");
    printf("%20s %9s %9s %9s %9s\n", "blocked DGEMM: t", "MFLOPS"
    		, "residuum", "res_alt", "sup-diff");
#endif


    hpc::util::WallTime<double> wallTime;

    for (Index i=MIN_SIZE; i<=MAX_SIZE; i+=granularity(i)) {
    	hpc::myblas::gecopy(i, i, C.data, C.incRow, C.incCol, C1.data, C1.incRow, C1.incCol);
    	hpc::myblas::gecopy(i, i, C.data, C.incRow, C.incCol, C2.data, C2.incRow, C2.incCol);
    	auto AV = A(0, 0, i, i);
    	auto BV = B(0, 0, i, i);
    	auto CV = C(0, 0, i, i);
    	auto C1V = C1(0, 0, i, i);
    	auto C2V = C2(0, 0, i, i);
        printf("%5ld ", i);

        wallTime.tic();
#if PRECISION == 1
    cblas_sgemm(
#elif PRECISION == 2
    cblas_dgemm(
#endif
		CblasColMajor, CblasNoTrans, CblasNoTrans,
	    i, i, i, alpha, AV.data, AV.incCol,
	    BV.data, BV.incCol,
	    beta, C1V.data, C1V.incCol);
        double t = wallTime.toc();
        printf("%20.4lf %9.2lf", t, 2.*i/1000*i/1000*i/t);

        wallTime.tic();
        hpc::myblas::gemm(i, i, i, alpha,
                      AV.data, AV.incRow, AV.incCol,
                      BV.data, BV.incRow, BV.incCol,
                      beta,
                      C2V.data, C2V.incRow, C2V.incCol);
        t = wallTime.toc();
        double res = hpc::matvec::residuum_gemm(alpha, AV, BV, beta, C1V, C2V);
        double res_alt = hpc::matvec::residuum_gemm_alt(alpha
        				, AV, BV, beta
        				, CV, C1V, C2V);
        double res_sup = hpc::matvec::max_rel_diff(C1V, C2V);
        printf("%20.4lf %9.2lf %9.1e %9.1e %9.1e", t, 2.*i/1000*i/1000*i/t, res, res_alt, res_sup);
        printf("\n");
    }
}
