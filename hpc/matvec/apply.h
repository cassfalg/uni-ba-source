#ifndef HPC_MATVEC_APPLY_H
#define HPC_MATVEC_APPLY_H 1

#include <type_traits>
#include <hpc/matvec/isgematrix.h>
#include <hpc/matvec/isdensevector.h>

namespace hpc { namespace matvec {

template <typename Vector, typename Func>
typename std::enable_if<IsDenseVector<Vector>::value,
         void>::type
apply(Vector& vector, Func func)
{
    using Index = typename Vector::Index;
    for (Index i = 0; i < vector.length; ++i) {
	func(vector(i), i);
    }
}

template <typename Vector, typename Func>
typename std::enable_if<IsDenseVector<Vector>::value,
         void>::type
apply(const Vector& vector, Func func)
{
    using Index = typename Vector::Index;
    for (Index i = 0; i < vector.length; ++i) {
	func(vector(i), i);
    }
}

template <typename MA, typename Func>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
apply(MA &A, Func func)
{
    typedef typename MA::Index    Index;

    if (A.incRow<A.incCol) {
        for (Index j=0; j<A.numCols; ++j) {
            for (Index i=0; i<A.numRows; ++i) {
                func(A(i,j), i, j);
            }
        }
    } else {
        for (Index i=0; i<A.numRows; ++i) {
            for (Index j=0; j<A.numCols; ++j) {
                func(A(i,j), i, j);
            }
        }
    }
}

template <typename MA, typename Func>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
apply(const MA &A, Func func)
{
    typedef typename MA::Index    Index;

    if (A.incRow<A.incCol) {
        for (Index j=0; j<A.numCols; ++j) {
            for (Index i=0; i<A.numRows; ++i) {
                func(A(i,j), i, j);
            }
        }
    } else {
        for (Index i=0; i<A.numRows; ++i) {
            for (Index j=0; j<A.numCols; ++j) {
                func(A(i,j), i, j);
            }
        }
    }
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_APPLY_H
