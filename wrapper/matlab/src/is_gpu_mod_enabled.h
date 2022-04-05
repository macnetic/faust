#ifndef IS_GPU_MOD_ENABLED
#define IS_GPU_MOD_ENABLED
#include "faust_constant.h"
template<FDevice DEV>
void _is_gpu_mod_enabled(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs);
#include "is_gpu_mod_enabled.hpp"
#endif
