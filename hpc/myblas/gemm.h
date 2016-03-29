#ifndef HPC_MYBLAS_GEMM_H
#define HPC_MYBLAS_GEMM_H 1

//#include <complex>
#include <type_traits>
#include <hpc/myblas/blocksize.h>
#include <hpc/myblas/gescal.h>
#include <hpc/myblas/geaxpy.h>
#include <hpc/myblas/kernels/ugemm_ref.h>
//#include <cstdlib>


namespace hpc { namespace myblas {

  template <typename Index, typename TA, typename T>
  void
  pack_A(Index mc, Index kc,
         const TA *A, Index incRowA, Index incColA,
         T *p)
  {
      Index MR = BlockSize<T>::MR;
      Index mp = (mc+MR-1) / MR;

      for (Index j=0; j<kc; ++j) {
          for (Index l=0; l<mp; ++l) {
              for (Index i0=0; i0<MR; ++i0) {
                  Index i  = l*MR + i0;
                  Index nu = l*MR*kc + j*MR + i0;
                  p[nu]   = (i<mc) ? A[i*incRowA+j*incColA]
                                   : T(0);
              }
          }
      }
  }

  template <typename Index, typename TB, typename T>
  void
  pack_B(Index kc, Index nc,
         const TB *B, Index incRowB, Index incColB,
         T *p)
  {
      Index NR = BlockSize<T>::NR;
      Index np = (nc+NR-1) / NR;

      for (Index l=0; l<np; ++l) {
          for (Index j0=0; j0<NR; ++j0) {
              for (Index i=0; i<kc; ++i) {
                  Index j  = l*NR+j0;
                  Index nu = l*NR*kc + i*NR + j0;
                  p[nu]   = (j<nc) ? B[i*incRowB+j*incColB]
                                   : T(0);
              }
          }
      }
  }

//-----------------------------------------------------------------------------

  template <typename Index, typename T, typename Beta, typename TC>
  void
  mgemm(Index mc, Index nc, Index kc,
        T alpha,
        const T *A, const T *B,
        Beta beta,
        TC *C, Index incRowC, Index incColC)
  {
      const Index MR = BlockSize<T>::MR;
      const Index NR = BlockSize<T>::NR;
      T C_[BlockSize<T>::MR*BlockSize<T>::NR];

      const Index mp  = (mc+MR-1) / MR;
      const Index np  = (nc+NR-1) / NR;
      const Index mr_ = mc % MR;
      const Index nr_ = nc % NR;

      for (Index j=0; j<np; ++j) {
          const Index nr = (j!=np-1 || nr_==0) ? NR : nr_;

          for (Index i=0; i<mp; ++i) {
              const Index mr = (i!=mp-1 || mr_==0) ? MR : mr_;

              if (mr==MR && nr==NR) {
                  ugemm(kc, alpha,
                        &A[i*kc*MR], &B[j*kc*NR],
                        beta,
                        &C[i*MR*incRowC+j*NR*incColC],
                        incRowC, incColC);
              } else {
                  std::fill_n(C_, MR*NR, T(0));
                  ugemm_0(kc, alpha,
                        &A[i*kc*MR], &B[j*kc*NR],
//                        T(0),
                        C_, Index(1), MR);
                  gescal(mr, nr, beta,
                         &C[i*MR*incRowC+j*NR*incColC],
                         incRowC, incColC);
                  geaxpy(mr, nr, T(1), C_, Index(1), MR,
                         &C[i*MR*incRowC+j*NR*incColC],
                         incRowC, incColC);
              }
          }
      }
  }

//-----------------------------------------------------------------------------

template <typename Index, typename Alpha,
         typename TA, typename TB,
         typename Beta,
         typename TC>
void
gemm(Index m, Index n, Index k,
     Alpha alpha,
     const TA *A, Index incRowA, Index incColA,
     const TB *B, Index incRowB, Index incColB,
     Beta beta,
     TC *C, Index incRowC, Index incColC)
{
    typedef typename std::common_type<Alpha, TA, TB>::type  T;

    const Index MC = BlockSize<T>::MC;
    const Index NC = BlockSize<T>::NC;
    const Index KC = BlockSize<T>::KC;

    const Index MR = BlockSize<T>::MR;
    const Index NR = BlockSize<T>::NR;

    const Index mb = (m+MC-1) / MC;
    const Index nb = (n+NC-1) / NC;
    const Index kb = (k+KC-1) / KC;

    const Index mc_ = m % MC;
    const Index nc_ = n % NC;
    const Index kc_ = k % KC;

    T *A_ = new T[MC*KC + MR];
    T *B_ = new T[KC*NC + NR];

    if (alpha==Alpha(0) || k==0) {
        gescal(m, n, beta, C, incRowC, incColC);
        return;
    }

    for (Index j=0; j<nb; ++j) {
        Index nc = (j!=nb-1 || nc_==0) ? NC : nc_;

        for (Index l=0; l<kb; ++l) {
            Index   kc  = (l!=kb-1 || kc_==0) ? KC : kc_;
            Beta beta_  = (l==0) ? beta : Beta(1);

            pack_B(kc, nc,
                   &B[l*KC*incRowB+j*NC*incColB],
                   incRowB, incColB,
                   B_);

            for (Index i=0; i<mb; ++i) {
                Index mc = (i!=mb-1 || mc_==0) ? MC : mc_;

                pack_A(mc, kc,
                       &A[i*MC*incRowA+l*KC*incColA],
                       incRowA, incColA,
                       A_);

                mgemm(mc, nc, kc,
                      T(alpha), A_, B_, beta_,
                      &C[i*MC*incRowC+j*NC*incColC],
                      incRowC, incColC);
            }
        }
    }
    delete[] A_;
    delete[] B_;
}

} } // namespace myblas, hpc

#endif // HPC_MYBLAS_GEMM_H
