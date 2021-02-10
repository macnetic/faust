#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_optimize_storage(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
			if(nrhs != 3)
				mexErrMsgTxt("optimize_storage mex error: invalid number of arguments.");
			bool timeCrit = (bool) mxGetScalar(prhs[2]);
			Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->optimize_storage(timeCrit);
			plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);

}

