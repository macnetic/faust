#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_disp(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 0 || nrhs != 2)
		mexErrMsgTxt("disp: Unexpected arguments");
	//core_ptr->display(); // doesn't work on windows,
	// matlab terminal doesn't receive the mex lib std output
	// we must use mexPrintf() to display content
	mexPrintf((core_ptr->to_string()+"\r\n").c_str());
}

