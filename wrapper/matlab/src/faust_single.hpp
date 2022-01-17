#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_single(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	Faust::TransformHelper<float,DEV>* th = core_ptr->template real<float>();
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<float,DEV> >(th);
}
