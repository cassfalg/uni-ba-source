#include <cstdio>
#include <hpc/util/walltime.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/random_init.h>
#include "blis/blis.h"

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

int main() {
    err_t blis_error = bli_init();
    if (blis_error != BLIS_SUCCESS) {
        printf("error initializing blis library\n");
    }
    using namespace hpc::matvec;
#if PRECISION == 1
    typedef float ElementType;
    printf("%7s ", "m=n=k");
    printf("%20s %9s\n", "bli_sgemm: t", "MFLOPS");
#elif PRECISION == 2
    typedef double ElementType;
    printf("%7s ", "m=n=k");
    printf("%20s %9s\n", "bli_dgemm: t", "MFLOPS");
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
        /* BLIS_NO_TRANSPOSE
           void bli_?gemm( trans_t transa,
           trans_t transb,
           dim_t   m,
           dim_t   n,
           dim_t   k,
           ctype*  alpha,
           ctype*  a, inc_t rsa, inc_t csa,
           ctype*  b, inc_t rsb, inc_t csb,
           ctype*  beta,
           ctype*  c, inc_t rsc, inc_t csc );*/
        #if PRECISION == 1
            bli_sgemm(
        #elif PRECISION == 2
            bli_dgemm(
        #endif
                BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
                i, i, i,
                &alpha,
                A.data, A.incCol, A.incRow,
                B.data, B.incCol, B.incRow,
                &beta,
                C.data, C.incCol, C.incRow);
        t = wallTime.toc();
        printf("%20.4lf %9.2lf\n", t, 2.*i/1000*i/1000*i/t);
    }
    blis_error = bli_finalize();
    if (blis_error != BLIS_SUCCESS) {
        printf("error finalizing blis library\n");
    }
}
