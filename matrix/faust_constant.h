#ifndef FAUST_CONSTANT_H
#define FAUST_CONSTANT_H

#ifdef FAUST_SINGLE
   typedef float faust_real;
#else
   typedef double faust_real;
#endif

const double FAUST_PRECISION = 0.0001;


#endif
