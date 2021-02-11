#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_optimize(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if(nrhs != 3)
		mexErrMsgTxt("optimize mex error: invalid number of arguments.");
	bool transp = (bool) mxGetScalar(prhs[2]);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->optimize(transp);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

template <typename SCALAR, FDevice DEV>
void set_FM_mul_mode(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	const int mode = (int) mxGetScalar(prhs[2]);
	core_ptr->set_FM_mul_mode(mode);
	return;
}

template <typename SCALAR, FDevice DEV>
void set_Fv_mul_mode(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	const int mode = (int) mxGetScalar(prhs[2]);
	core_ptr->set_Fv_mul_mode(mode);
	return;
}
