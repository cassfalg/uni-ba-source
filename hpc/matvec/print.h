#ifndef HPC_MATVEC_PRINT_H
#define HPC_MATVEC_PRINT_H 1

#include <complex>
#include <type_traits>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename T>
struct Format
{
};

template <>
struct Format<float>
{
    static const char *
    value()
    {
        return " %10.1f";
    }
};

template <>
struct Format<double>
{
    static const char *
    value()
    {
        return " %10.1lf";
    }
};

template <>
struct Format<std::complex<double> >
{
    static const char *
    value()
    {
        return " (%10.1lf, %10.1lf)";
    }
};

//------------------------------------------------------------------------------

template <typename T>
void
print_value(T value) {
   std::printf(Format<T>::value(), value);
}

template <typename T>
void
print_value(const std::complex<T> &value) {
   std::printf(Format<std::complex<T> >::value(), value.real(), value.imag());
}

template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
print(const MA &A, const char *name = "")
{
    typedef typename MA::Index    Index;

    if (*name) {
        printf("%s = \n", name);
    }
    for (Index i=0; i<A.numRows; ++i) {
        for (Index j=0; j<A.numCols; ++j) {
            print_value(A(i,j));
        }
        printf("\n");
    }
    printf("\n");
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_PRINT_H
