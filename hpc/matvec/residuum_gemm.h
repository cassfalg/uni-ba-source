#ifndef HPC_MATVEC_RESIDUUM_GEMM_H
#define HPC_MATVEC_RESIDUUM_GEMM_H 1

#include <random>
#include <complex>
#include <limits>
#include <cmath>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

  //
  //  Compute 1-norm of difference A - B
  //
template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value, typename MA::ElementType>::type
asumDiff(const MA &A, const MA &B) {
	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	ElementType diff = ElementType(0);

	for (Index i=0; i<A.numRows; ++i) {
		for (Index j=0; j<A.numCols; ++j) {
			diff += std::abs(A(i,j) - B(i,j));
		}
	}
	return diff;
}


  //
  //  Compute sup-norm of difference A - B
  //
template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value, typename MA::ElementType>::type
asupDiff(const MA &A, const MA &B)
{
	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	ElementType maxDiff = ElementType(0);
	ElementType curDiff = ElementType(0);

	for (Index i=0; i<A.numRows; ++i) {
		for (Index j=0; j<A.numCols; ++j) {
			curDiff = std::abs(A(i, j) - B(i, j));
			maxDiff = (curDiff>maxDiff)?curDiff:maxDiff;
		}
	}
	return maxDiff;
}

  //
  // Compute 1-norm
  //
template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value, typename MA::ElementType>::type
asum(const MA &A) {
	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	ElementType result(0);

	for (Index i=0; i<A.numRows; ++i) {
		for (Index j=0; j<A.numCols; ++j) {
			result += std::abs(A(i, j));
		}
	}
	return result;
}

  //
  // Compute sup-norm
  //
template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value, typename MA::ElementType>::type
asup(const MA &A) {
	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	ElementType max(0);
	ElementType cur(0);
	for (Index i=0; i<A.numRows; ++i) {
		for (Index j=0; j<A.numCols; ++j) {
			cur = std::abs(A(i, j));
			max = (cur > max)?cur:max;
		}
	}
	return max;
}

  //
  //  Compute residual for matrix-product with 1-Norm
  //
template <typename MA>
typename MA::ElementType
residuum_gemm(typename MA::Index n, typename MA::ElementType alpha
		, typename MA::ElementType aNorm, typename MA::ElementType bNorm
		, typename MA::ElementType beta
		, typename MA::ElementType cNorm, MA &C1, MA &C2)
{
	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	Index maxDim = std::max(std::max(C1.numRows, n), C1.numCols);
	ElementType aDiff = asumDiff(C1, C2);
      // Using eps for double gives upper bound in case elements have lower
      // precision.
	ElementType eps = std::numeric_limits<ElementType>::epsilon();
	ElementType res = aDiff/(aNorm*std::abs(alpha)
						*bNorm
						*cNorm*std::max(ElementType(1),std::abs(beta))
						*eps*maxDim);
	return res;
}

  //
  //  Compute residual for matrix-product with sup norm
  // GEMM: C = beta * C + alpha * A * B
template <typename MA>
typename MA::ElementType
max_rel_diff(const MA &A, const MA &B)
{

	typedef typename MA::Index Index;
	typedef typename MA::ElementType ElementType;
	ElementType size(0);
	ElementType diff(0);
	ElementType max(0);
	ElementType cur(0);
	ElementType tresh(std::numeric_limits<ElementType>::epsilon()*ElementType(100));
	for (Index i=0; i<A.numRows; ++i) {
		for (Index j=0; j<A.numCols; ++j) {
			diff = std::abs(A(i, j) - B(i, j));
			size = std::abs(A(i, j) + B(i, j))/2;
			cur = (size > tresh)?diff/size:diff;
			max = (cur > max)?cur:max;
		}
	}
	return max;}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_RESIDUUM_GEMM_H
