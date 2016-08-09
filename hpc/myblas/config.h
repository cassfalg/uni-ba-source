#ifndef HPC_MYBLAS_CONFIG_H
#define HPC_MYBLAS_CONFIG_H 1


//
// default
//
#ifndef DEFAULT_BLOCKSIZE_MR
#define DEFAULT_BLOCKSIZE_MR 4
#endif

#ifndef DEFAULT_BLOCKSIZE_NR
#define DEFAULT_BLOCKSIZE_NR 4
#endif

#ifndef DEFAULT_BLOCKSIZE_MC
#define DEFAULT_BLOCKSIZE_MC 64
#endif

#ifndef DEFAULT_BLOCKSIZE_KC
#define DEFAULT_BLOCKSIZE_KC 64
#endif

#ifndef DEFAULT_BLOCKSIZE_NC
#define DEFAULT_BLOCKSIZE_NC 256
#endif

//
// single precision(float)
//
#ifndef S_BLOCKSIZE_MR
#define S_BLOCKSIZE_MR 4
#endif

#ifndef S_BLOCKSIZE_NR
#define S_BLOCKSIZE_NR 4
#endif

#ifndef S_BLOCKSIZE_MC
#define S_BLOCKSIZE_MC 320
#endif

#ifndef S_BLOCKSIZE_KC
#define S_BLOCKSIZE_KC 320
#endif

#ifndef S_BLOCKSIZE_NC
#define S_BLOCKSIZE_NC 5120
#endif

//
// double precision(double)
//
#ifndef D_BLOCKSIZE_MR
#define D_BLOCKSIZE_MR 4
#endif

#ifndef D_BLOCKSIZE_NR
#define D_BLOCKSIZE_NR 4
#endif

#ifndef D_BLOCKSIZE_MC
#define D_BLOCKSIZE_MC 192
#endif

#ifndef D_BLOCKSIZE_KC
#define D_BLOCKSIZE_KC 192
#endif

#ifndef D_BLOCKSIZE_NC
#define D_BLOCKSIZE_NC 3072
#endif

//
// complex single precision(std::complex<float>)
//
#ifndef C_BLOCKSIZE_MR
#define C_BLOCKSIZE_MR 4
#endif

#ifndef C_BLOCKSIZE_NR
#define C_BLOCKSIZE_NR 8
#endif

#ifndef C_BLOCKSIZE_MC
#define C_BLOCKSIZE_MC 256
#endif

#ifndef C_BLOCKSIZE_KC
#define C_BLOCKSIZE_KC 256
#endif

#ifndef C_BLOCKSIZE_NC
#define C_BLOCKSIZE_NC 4096
#endif

//
// complex double precision(std::complex<double>)
//
#ifndef Z_BLOCKSIZE_MR
#define Z_BLOCKSIZE_MR 4
#endif

#ifndef Z_BLOCKSIZE_NR
#define Z_BLOCKSIZE_NR 4
#endif

#ifndef Z_BLOCKSIZE_MC
#define Z_BLOCKSIZE_MC 256
#endif

#ifndef Z_BLOCKSIZE_KC
#define Z_BLOCKSIZE_KC 128
#endif

#ifndef Z_BLOCKSIZE_NC
#define Z_BLOCKSIZE_NC 4096
#endif

#endif // HPC_MYBLAS_BLOCKSIZE_H
