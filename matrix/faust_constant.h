#ifndef FAUST_CONSTANT_H
#define FAUST_CONSTANT_H

#ifdef FAUST_SINGLE
   typedef float faust_real;
#else
   typedef double faust_real;
#endif

typedef unsigned long int faust_unsigned_int;
typedef long int faust_int;

const double FAUST_PRECISION = 0.0001;
const faust_unsigned_int nbr_iter_norm = 100; 


#endif
