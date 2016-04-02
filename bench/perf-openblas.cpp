#include <cstdio>
#include <hpc/util/walltime.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/random_init.h>
#include <cblas.h>

#ifndef MIN_SIZE
#define MIN_SIZE 200
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 4000
#endif
#ifndef GRANULARITY
#define GRANULARITY 200u
#endif

int main() {
  openblas_set_num_threads(1);
  using namespace hpc::matvec;
  typedef double ElementType;
  typedef int Index; //openblas uses int
  GeMatrix<ElementType, Index> A = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
  GeMatrix<ElementType, Index> B = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
  GeMatrix<ElementType, Index> C = GeMatrix<ElementType, Index>(MAX_SIZE, MAX_SIZE);
  randomInit(A);
  randomInit(B);
  randomInit(C);
  const ElementType alpha(1.5);
  const ElementType beta(2.5);

  printf("%7s ", "m=n=k");
  printf("%20s %9s", "OpenBLAS: t", "MFLOPS\n");
  hpc::util::WallTime<double> wallTime;
  double t;
  for (Index i=MIN_SIZE; i<=MAX_SIZE; i+=GRANULARITY) {
      printf("  %5d ", i);
	wallTime.tic();
	/*void cblas_dgemm(
		   CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
		   const int M, const int N, const int K,
		   const double alpha, const double *A, const int lda,
		   const double *B, const int ldb,
		   const double beta, double *C, const int ldc);*/
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		    i, i, i, alpha, &A(0,0), A.incCol,
		    B.data, B.incCol,
		    beta, C.data, C.incCol);
	t = wallTime.toc();
      printf("%20.4lf %9.2lf\n", t, 2.*i/1000*i/1000*i/t);
  }
}
