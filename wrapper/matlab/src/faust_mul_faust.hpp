#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_mul_faust(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs != 3)
		mexErrMsgTxt("mul_faust: Unexpected arguments");
	Faust::TransformHelper<SCALAR,DEV>* right_term = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[2]);
	if(right_term->getNbRow() != core_ptr->getNbCol()) mexErrMsgTxt("The dimensions of the two Fausts must agree.");
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->multiply(right_term);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
#ifdef FAUST_VERBOSE
	printf("new ptr for Faust multiply's TransformHelper obj=%p\n", th);
	printf("multiply(\"faust\")\n");
	printf("isConjugate()=%d\n",th->isConjugate());
	printf("isTransposed()=%d\n", th->isTransposed());
#endif
}

