#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV, typename FPP>
void faust_mul_scalar(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs != 3)
	{
		mexErrMsgTxt("mul_scalar: Unexpected arguments");
	}
	//complex Faust times real scalar OK
	//complex Faust times complex oK
	//real Faust times complex NOT OK
	//real Faust times real OK
	if(mxIsComplex(prhs[2]) && (typeid(SCALAR) == typeid(double) || typeid(SCALAR) == typeid(float)))
		mexErrMsgTxt("Multiplying a real Faust by a complex scalar is not yet implemented.");

	SCALAR scalar;
	mxArray2Scalar<FPP>(prhs[2], &scalar);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->multiply(scalar);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

