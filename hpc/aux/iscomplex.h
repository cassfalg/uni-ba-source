#ifndef HPC_AUX_ISCOMPLEX_H
#define HPC_AUX_ISCOMPLEX_H 1

#include <complex>
#include <type_traits>

namespace hpc { namespace aux {

template <typename Any>
struct IsComplex
{
    struct Two {
        char a, b;
    };

    template <typename T>
        static char
        checker(const std::complex<T> &);

    template <typename T>
        static Two
        checker(const T &);

    static Any        any;
    static const bool value = sizeof(checker(any)) == 1;
};



} } // namespace aux, hpc

#endif // HPC_AUX_ISCOMPLEX_H
