#ifndef FAUST_CONSTANT_H
#define FAUST_CONSTANT_H

///if defined, compilation with openblas for dense algebra
#define __GEMM_WITH_OPENBLAS__

/*! \brief If defined, faust uses single precision for the the algebra (Faust::MatDense, Faust::MatSparse, Faust::Transform). <br>
* If not, it's used double precision.
*\warning Only useful for the mexfunction (MATLAB wrapper) since the Faust::MatDense ... classes are template
*/
/* #undef FAUST_SINGLE */

//if defined, print debug message (useful for debugging)
/* #undef FAUST_VERBOSE */

//if defined, print timers (profile the code)
/* #undef __COMPILE_TIMERS__ */


/// Value of FFPP (faust floating point precision)
#ifdef FAUST_SINGLE
	typedef float FFPP; ///single precision
#else
	typedef double FFPP; ///double precision
#endif

enum Device
{
   Cpu, /* Central Process Unit */
#ifdef __COMPILE_GPU__
   Gpu, /* Graphic Process Unit */
#endif
};
typedef unsigned long int faust_unsigned_int;
typedef long int faust_int;

const double FAUST_PRECISION = 0.0001;
const faust_unsigned_int nbr_iter_norm = 100;

///// #define FAUST_GPU_ALREADY_SET -3
///// #define FAUST_DEFAULT_CUDA_DEVICE 0

#endif
