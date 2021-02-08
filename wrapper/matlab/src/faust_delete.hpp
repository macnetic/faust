#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_delete(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	// Destroy the C++ object
	destroyObject<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	// Warn if other commands were ignored
	if (nlhs != 0 || nrhs != 2)
		mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
}

