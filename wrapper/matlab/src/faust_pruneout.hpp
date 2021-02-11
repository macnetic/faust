#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_pruneout(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	int nnz_tres = 0;
	int npasses = -1;
	bool only_forward = false;
	if (nlhs == 1 || nrhs == 5)
	{
		nnz_tres = (int) mxGetScalar(prhs[2]);
		npasses = (int) mxGetScalar(prhs[3]);
		only_forward = (bool) mxGetScalar(prhs[4]);
	}
	else
		mexErrMsgTxt("pruneout mex error: number of arguments must be 3.");
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->pruneout(nnz_tres, npasses, only_forward);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

