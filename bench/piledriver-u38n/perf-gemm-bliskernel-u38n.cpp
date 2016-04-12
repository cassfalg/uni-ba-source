#include <cstdio>
#include <hpc/util/walltime.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/random_init.h>

#ifndef MIN_SIZE
#define MIN_SIZE 200
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 4000
#endif
#ifndef GRANULARITY
#define GRANULARITY 200u
#endif

//    1: single precision
//    2: double precision
#ifndef PRECISION
#define PRECISION 2
#endif

#include "mygemm-blisugem-u38n.h"

int main() {
/*    err_t blis_error = bli_init();
    if (blis_error != BLIS_SUCCESS) {
        printf("error initializing blis library\n");
    }
*/
    using namespace hpc::matvec;
#if PRECISION == 1
    typedef float ElementType;
    printf("%7s ", "m=n=k");
    printf("%20s %9s\n", "hpc_sgemm:t", "MFLOPS");
#elif PRECISION == 2
    typedef double ElementType;
    printf("%7s ", "m=n=k");
    printf("%20s %9s\n", "hpc_dgemm:t", "MFLOPS");
#endif
    typedef gint_t Index; //openblas uses int
    GeMatrix<ElementType, Index> A = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
    GeMatrix<ElementType, Index> B = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
    GeMatrix<ElementType, Index> C = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
    randomInit(A);
    randomInit(B);
    randomInit(C);
    ElementType alpha(1.5);
    ElementType beta(2.5);

    hpc::util::WallTime<double> wallTime;
    double t;
    for (Index i=MIN_SIZE; i<=MAX_SIZE; i+=GRANULARITY) {
        printf("  %5ld ", i);
        wallTime.tic();
    hpc::myblas::gemm(i, i, i, alpha,
		       A.data, A.incRow, A.incCol,
		       B.data, B.incRow, B.incCol,
		       beta,
		       C.data, C.incRow, C.incCol);
        t = wallTime.toc();
        printf("%20.4lf %9.2lf\n", t, 2.*i/1000*i/1000*i/t);
    }
/*    blis_error = bli_finalize();
    if (blis_error != BLIS_SUCCESS) {
        printf("error finalizing blis library\n");
    }
*/
}
