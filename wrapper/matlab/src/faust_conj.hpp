#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_conj(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->conjugate();
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
#ifdef FAUST_VERBOSE
	printf("new ptr for Faust conjugate's TransformHelper obj=%p\n", th);
	printf("conjugate()\n");
	printf("isConjugate()=%d\n",th->isConjugate());
#endif

}

