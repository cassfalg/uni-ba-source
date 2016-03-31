#include <complex>
#include <cstdio>
#include <hpc/myblas/gecopy.h>
#include <hpc/myblas/gemm.h>
#include <hpc/refblas/gemm.h>
#include <hpc/matvec/residuum_gemm.h>
#include <hpc/matvec/random_init.h>
#include <hpc/util/walltime.h>

int
main()
{
    // Element type for A, B and C
    typedef float                TA;
    typedef float                TB;
//    typedef std::complex<double> TC;
    typedef float                TC;

    // Type for scalar factors
    typedef float   Alpha;
    typedef float   Beta;

    // Max dimensions of A, B and C 
    long max_m  = 7000;
    long max_n  = 7000;
    long max_k  = 7000;

    // Storage order of A
    long incRowA = 1;
    long incColA = max_m;

    // Storage order of B
    long incRowB = max_n;
    long incColB = 1;

    // Storage order of C
    long incRowC = max_n;
    long incColC = 1;

    // Allocate A, B, C1, C2
    TA *A  = new TA[max_m*max_k];
    TB *B  = new TB[max_k*max_n];
    TC *C1 = new TC[max_m*max_n];
    TC *C2 = new TC[max_m*max_n];

    // Init all matrices
    hpc::matvec::randomInit(max_m, max_k, A, incRowA, incColA);
    hpc::matvec::randomInit(max_k, max_n, B, incRowB, incColB);
    hpc::matvec::randomInit(max_m, max_n, C1, incRowC, incColC);
    hpc::myblas::gecopy(max_m, max_n,
                         C1, incRowC, incColC,
                         C2, incRowC, incColC);

    // Init scalar factors
    const Alpha alpha(1.5);
    const Beta  beta(2.5);

    // Header for benchmark
    printf("%5s %5s %5s ", "m", "n", "k");
    printf("%20s %9s", "refColMajor: t", "MFLOPS");
    printf("%20s %9s %9s", "blocked GEMM: t", "MFLOPS", "diff");
    printf("\n");

    hpc::util::WallTime<double> wallTime;

    for (long m=200, n=200, k=200;
         m <=max_m && n<=max_n && k<=max_k;
         m+=100, n+=100, k+=100)
    {
        printf("%5ld %5ld %5ld ", m, n, k);

        wallTime.tic();
        hpc::refblas::gemm(m, n, k, alpha,
                      A, incRowA, incColA,
                      B, incRowB, incColB,
                      beta,
                      C1, incRowC, incColC);
        double t = wallTime.toc();
        printf("%20.4lf %9.2lf", t, 2.*m/1000*n/1000*k/t);

        wallTime.tic();
        hpc::myblas::gemm(m, n, k, alpha,
                           A, incRowA, incColA,
                           B, incRowB, incColB,
                           beta,
                           C2, incRowC, incColC);
        t = wallTime.toc();
        double res = hpc::matvec::residuum_gemm(m, n, k, alpha,
                              A, incRowA, incColA,
                              B, incRowB, incColB,
                              beta,
                              C1, incRowC, incColC,
                              C2, incRowC, incColC);

        printf("%20.4lf %9.2lf %9.1e", t, 2.*m/1000*n/1000*k/t, res);
        printf("\n");
    }
}
