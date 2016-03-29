#ifndef HPC_MYBLAS_BLOCKSIZE_H
#define HPC_MYBLAS_BLOCKSIZE_H 1

#include <complex>
#include <hpc/myblas/config.h>

namespace hpc { namespace myblas {

template <typename T>
struct BlockSize
{
    static const int MC = DEFAULT_BLOCKSIZE_MC;
    static const int KC = DEFAULT_BLOCKSIZE_KC;
    static const int NC = DEFAULT_BLOCKSIZE_NC;
    static const int MR = DEFAULT_BLOCKSIZE_MR;
    static const int NR = DEFAULT_BLOCKSIZE_NR;


    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<float>
{
    static const int MC = S_BLOCKSIZE_MC;
    static const int KC = S_BLOCKSIZE_KC;
    static const int NC = S_BLOCKSIZE_NC;
    static const int MR = S_BLOCKSIZE_MR;
    static const int NR = S_BLOCKSIZE_NR;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<double>
{
    static const int MC = D_BLOCKSIZE_MC;
    static const int KC = D_BLOCKSIZE_KC;
    static const int NC = D_BLOCKSIZE_NC;
    static const int MR = D_BLOCKSIZE_MR;
    static const int NR = D_BLOCKSIZE_NR;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<std::complex<float> >
{
    static const int MC = C_BLOCKSIZE_MC;
    static const int KC = C_BLOCKSIZE_KC;
    static const int NC = C_BLOCKSIZE_NC;
    static const int MR = C_BLOCKSIZE_MR;
    static const int NR = C_BLOCKSIZE_NR;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

template <>
struct BlockSize<std::complex<double> >
{
    static const int MC = Z_BLOCKSIZE_MC;
    static const int KC = Z_BLOCKSIZE_KC;
    static const int NC = Z_BLOCKSIZE_NC;
    static const int MR = Z_BLOCKSIZE_MR;
    static const int NR = Z_BLOCKSIZE_NR;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};

} } // namespace myblas, hpc

#endif // HPC_MYBLAS_BLOCKSIZE_H
