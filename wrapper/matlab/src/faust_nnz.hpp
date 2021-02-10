#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_nnz(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
		// return the number of nonzeros coefficients of the matrix
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs > 1 || nrhs != 2)
	{
		mexErrMsgTxt("nnz : incorrect number of arguments.");
	}

	long long int nnz = core_ptr->get_total_nnz();
	plhs[0]=mxCreateDoubleScalar(nnz);
}

