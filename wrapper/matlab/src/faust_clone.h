#ifndef __MEX_FAUST_CLONE__
#define __MEX_FAUST_CLONE__
template <typename SCALAR, FDevice DEV>
void faust_clone_cpu2gpu(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs);
template <typename SCALAR, FDevice DEV>
void faust_clone_gpu2cpu(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs);
#include "faust_clone.hpp"
#endif
