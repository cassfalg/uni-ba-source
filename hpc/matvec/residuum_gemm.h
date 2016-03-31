#ifndef HPC_MATVEC_RESIDUUM_GEMM_H
#define HPC_MATVEC_RESIDUUM_GEMM_H 1

#include <random>
#include <complex>
#include <limits>
#include <cmath>
#include <hpc/matvec/gematrix.h>

namespace hpc { namespace matvec {

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
  residuum_gemm(Index m, Index n, Index k,
           Alpha alpha,
           const TA *A, Index incRowA, Index incColA,
           const TB *B, Index incRowB, Index incColB,
           Beta beta,
           TC *C1, Index incRowC1, Index incColC1,
           TC *C2, Index incRowC2, Index incColC2)
  {
      double aNorm = asum(m, k, A, incRowA, incColA) * std::abs(alpha);
      double bNorm = asum(k, n, B, incRowB, incColB);
      double cNorm = asum(m, n, C2, incRowC2, incColC2) * std::max(Beta(1),std::abs(beta));
      double aDiff = asumDiff(m, n,
                              C1, incRowC1, incColC1,
                              C2, incRowC2, incColC2);
      // Using eps for double gives upper bound in case elements have lower
      // precision.
      double eps = std::numeric_limits<double>::epsilon();
      double res = aDiff/(aNorm*bNorm*cNorm*eps*std::max(std::max(m,n),k));
      return res;
  }

} } // namespace matvec, hpc

#endif // HPC_MATVEC_RESIDUUM_GEMM_H
