#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_copy(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs > 1)
		mexErrMsgTxt("copy : too many output arguments");
	else
	{
		TransformHelper<SCALAR,DEV>* th = core_ptr->clone();
		plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
	}

}

