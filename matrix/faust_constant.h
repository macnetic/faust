#ifndef FAUST_CONSTANT_H
#define FAUST_CONSTANT_H

// #ifdef FAUST_SINGLE
   // typedef float faust_real;
// #else
   // typedef double faust_real;
// #endif

#ifdef __FAUST_SINGLE
   typedef float faust_real;
#else
   typedef double faust_real;
#endif

#ifdef __FAUST_SINGLE__
	#define IS_SINGLE_DEFINED 1
#else
	#define IS_SINGLE_DEFINED 0
#endif



typedef unsigned long int faust_unsigned_int;
typedef long int faust_int;

const double FAUST_PRECISION = 0.0001;
const faust_unsigned_int nbr_iter_norm = 100; 

#define FAUST_GPU_ALREADY_SET -3
#define FAUST_DEFAULT_CUDA_DEVICE 0

#endif
