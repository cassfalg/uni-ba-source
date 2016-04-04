/*
 * perf-blis.cpp
 *
 *  Created on: 26.03.2016
 *      Author: chris
 *  Test gemm performance of blis
 */

#include <iomanip>
#include <iostream>

#include "utl_walltime.h"

#include "blis.h"


#ifndef MIN_SIZE
#define MIN_SIZE 100
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 800
#endif
#ifndef GRANULARITY
#define GRANULARITY 100u
#endif

/* BLIS_NO_TRANSPOSE
 * void bli_?gemm( trans_t transa,
                trans_t transb,
                dim_t   m,
                dim_t   n,
                dim_t   k,
                ctype*  alpha,
                ctype*  a, inc_t rsa, inc_t csa,
                ctype*  b, inc_t rsb, inc_t csb,
                ctype*  beta,
                ctype*  c, inc_t rsc, inc_t csc );*/

int
main ()
{
  err_t bli_init( void );
//  double* A = new double[MAX_SIZE * MAX_SIZE];
//  double* B = new double[MAX_SIZE * MAX_SIZE];
//  double* C = new double[MAX_SIZE * MAX_SIZE];

  double *A = (double *) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));
  double *B = (double *) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));
  double *C = (double *) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));
  double alpha = 1.0; double beta = 0.0;

//  float *A = (float *) malloc(MAX_SIZE * MAX_SIZE * sizeof(float));
//  float *B = (float *) malloc(MAX_SIZE * MAX_SIZE * sizeof(float));
//  float *C = (float *) malloc(MAX_SIZE * MAX_SIZE * sizeof(float));
//  float alpha = 1.0; float beta = 0.0;

    std::cout << "matrix_dim     time in ns" << std::endl;
  double start = walltime();
  for (long int m = MIN_SIZE; m < MAX_SIZE; m += GRANULARITY)
    {
      bli_dgemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
		m, m, m,
		&alpha,
		A, 1, m,
		B, 1, m,
		&beta,
		C, 1, m);
      std::cout << " " << std::setw (5) << m;
      std::cout << "       " << std::fixed << std::setw (10) << std::setprecision (5)
	  << walltime() - start;
      std::cout << std::endl;
      std::cout.flush ();
    }
  err_t bli_finalize( void );
}
