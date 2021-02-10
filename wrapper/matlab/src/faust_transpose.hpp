#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_transpose(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->transpose();
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
#ifdef FAUST_VERBOSE
	printf("new ptr for Faust transpose's TransformHelper obj=%p\n", th);
	printf("transpose()\n");
	printf("isTransposed()=%d\n",th->isTransposed());
#endif

}

