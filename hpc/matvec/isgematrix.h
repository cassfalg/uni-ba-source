#ifndef HPC_MATVEC_ISGEMATRIX_H
#define HPC_MATVEC_ISGEMATRIX_H 1

#include <cassert>
#include <type_traits>
#include <hpc/matvec/gematrix.h>

namespace hpc { namespace matvec {

template <typename Any>
struct IsGeMatrix_
{
    static constexpr bool value = false;
};

template <typename T, typename I>
struct IsGeMatrix_<GeMatrix<T,I> >
{
    static constexpr bool value = true;
};

template <typename T, typename I>
struct IsGeMatrix_<GeMatrixView<T,I> >
{
    static constexpr bool value = true;
};

template <typename T, typename I>
struct IsGeMatrix_<GeMatrixConstView<T,I> >
{
    static constexpr bool value = true;
};

template <typename Any_>
struct IsGeMatrix
{
    typedef typename std::remove_reference<Any_>::type  Any;
    static constexpr bool value = IsGeMatrix_<Any>::value;
};

} } // namespace matvec, hpc

#endif // HPC_MATVEC_ISGEMATRIX_H
