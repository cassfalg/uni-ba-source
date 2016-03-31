#ifndef HPC_MATVEC_RANDOM_INIT_H
#define HPC_MATVEC_RANDOM_INIT_H 1

#include <random>
#include <complex>
#include <hpc/matvec/gematrix.h>

namespace hpc { namespace matvec {


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


} } // namespace matvec, hpc

#endif // HPC_MATVEC_RANDOM_INIT_H
