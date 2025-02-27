#include "class_handle.hpp"
#include "faust_TransformHelper.h"
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#endif

template <typename SCALAR, FDevice DEV>
void faust_opt_butterfly(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	auto F = convertMat2Ptr<Faust::TransformHelper<SCALAR, DEV>>(prhs[1]);
	auto oF = Faust::TransformHelper<SCALAR, DEV>::optButterflyFaust(F);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV>>(oF);
}

