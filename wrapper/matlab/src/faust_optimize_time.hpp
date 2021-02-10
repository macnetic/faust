#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_optimize_time(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if(nrhs != 5)
		mexErrMsgTxt("optimize_time mex error: invalid number of arguments.");
	if(nlhs != 1 && nlhs != 0)
		mexErrMsgTxt("optimize_time mex error: this function doesn't return more than one argument.");
	bool transp = (bool) mxGetScalar(prhs[2]);
	bool inplace = (bool) mxGetScalar(prhs[3]);
	int nsamples = (int) mxGetScalar(prhs[4]);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->optimize_time(transp, inplace, nsamples);
	if(inplace /*th == nullptr*/)
		return;
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

