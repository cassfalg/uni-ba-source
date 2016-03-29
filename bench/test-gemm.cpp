#include <chrono>
#include <complex>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <hpc/myblas/gecopy.h>
#include <hpc/myblas/gemm.h>
#include <hpc/refblas/gemm.h>


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


//
//  Compute 1-norm of difference A - B
//
template <typename Index, typename T>
double
asumDiff(Index m, Index n,
         const T *A, Index incRowA, Index incColA,
         const T *B, Index incRowB, Index incColB)
{
    double diff = 0;

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            diff += std::abs(A[i*incRowA+j*incColA] - B[i*incRowB+j*incColB]);
        }
    }
    return diff;
}

//
// Compute 1-norm
//
template <typename Index, typename T>
double
asum(Index m, Index n,
     const T *A, Index incRowA, Index incColA)
{
    double result = 0;

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            result += std::abs(A[i*incRowA+j*incColA]);
        }
    }
    return result;
}

//
//  Compute residual for matrix-product
//
template <typename Index, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
double
res_gemm(Index m, Index n, Index k,
         Alpha alpha,
         const TA *A, Index incRowA, Index incColA,
         const TB *B, Index incRowB, Index incColB,
         Beta beta,
         TC *C1, Index incRowC1, Index incColC1,
         TC *C2, Index incRowC2, Index incColC2)
{
    double aNorm = asum(m, k, A, incRowA, incColA) * std::abs(alpha);
    double bNorm = asum(k, n, B, incRowB, incColB);
    double cNorm = asum(m, n, C2, incRowC2, incColC2);
    double aDiff = asumDiff(m, n,
                            C1, incRowC1, incColC1,
                            C2, incRowC2, incColC2);
    // Using eps for double gives upper bound in case elements have lower
    // precision.
    double eps = std::numeric_limits<double>::epsilon();
    double res = aDiff/(aNorm*bNorm*cNorm*eps*std::max(std::max(m,n),k));
    return res;
}

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
    randomInit(max_m, max_k, A, incRowA, incColA);
    randomInit(max_k, max_n, B, incRowB, incColB);
    randomInit(max_m, max_n, C1, incRowC, incColC);
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

    WallTime<double> wallTime;

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
        double res = res_gemm(m, n, k, alpha,
                              A, incRowA, incColA,
                              B, incRowB, incColB,
                              beta,
                              C1, incRowC, incColC,
                              C2, incRowC, incColC);

        printf("%20.4lf %9.2lf %9.1e", t, 2.*m/1000*n/1000*k/t, res);
        printf("\n");
    }
}
