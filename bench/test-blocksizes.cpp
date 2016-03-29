#include <chrono>
#include <complex>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <hpc/myblas/gemm.h>
#include <hpc/refblas/gemm.h>
#include <hpc/myblas/blocksize.h>


//
//  Random initializer for general matrices: real and complex valued
//
template <typename Index, typename T>
void
randomInit(Index m, Index n, T *A, Index incRowA, Index incColA)
{
    std::random_device                  random;
    std::default_random_engine          mt(random());
    std::uniform_real_distribution<T>   uniform(-100,100);

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            A[i*incRowA+j*incColA] = uniform(mt);
        }
    }
}

template <typename Index, typename T>
void
randomInit(Index m, Index n, std::complex<T> *A, Index incRowA, Index incColA)
{
    std::random_device                  random;
    std::default_random_engine          mt(random());
    std::uniform_real_distribution<T>   uniform(-100,100);

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            A[i*incRowA+j*incColA] = std::complex<T>(uniform(mt), uniform(mt));
        }
    }
}

//
// Timer
//
template <typename T>
struct WallTime
{
    void
    tic()
    {
        t0 = std::chrono::high_resolution_clock::now();
    }

    T
    toc()
    {
        using namespace std::chrono;

        elapsed = high_resolution_clock::now() - t0;
        return duration<T,seconds::period>(elapsed).count();
    }

    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::duration   elapsed;
};


int
main()
{
    // Element type for A, B and C
    typedef double                TA;
    typedef double                TB;
//    typedef std::complex<double> TC;
    typedef double                TC;

    // Type for scalar factors
    typedef double   Alpha;
    typedef double   Beta;

    // Max dimensions of A, B and C 
    const long max_m  = 1200;
    const long max_n  = 1200;
    const long max_k  = 1200;

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

    // Init all matrices
    randomInit(max_m, max_k, A, incRowA, incColA);
    randomInit(max_k, max_n, B, incRowB, incColB);
    randomInit(max_m, max_n, C1, incRowC, incColC);

    // Init scalar factors
    const Alpha alpha(1.5);
    const Beta  beta(2.5);

    // Header for benchmark
//    printf("%5s %5s %5s %5s %5s", "MC", "KC", "NC", "MR", "NR");
//    printf("%20s %9s", "blocked GEMM: t", "MFLOPS");
//    printf("\n");

    WallTime<double> wallTime;


    printf("%5d %5d %5d %5d %5d ", D_BLOCKSIZE_MC, D_BLOCKSIZE_KC,
	   D_BLOCKSIZE_NC, D_BLOCKSIZE_MR, D_BLOCKSIZE_NR);

    wallTime.tic();
    hpc::myblas::gemm(max_m, max_n, max_k, alpha,
		       A, incRowA, incColA,
		       B, incRowB, incColB,
		       beta,
		       C1, incRowC, incColC);
    double t = wallTime.toc();

    printf("%20.4lf %9.2lf ", t, 2.*max_m/1000*max_n/1000*max_k/t);
    printf("\n");

}
