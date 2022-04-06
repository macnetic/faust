#include "faust2Mx.h"
#include "class_handle.hpp"
#include "faust_TransformHelper.h"


template <typename SCALAR, FDevice DEV>
void faust_is_all_sparse(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs != 1)
	{
		mexErrMsgTxt("is_all_sparse: incorrect number of output arguments.");
	}
	if (nrhs != 4)
	{
		mexErrMsgTxt("is_all_sparse: incorrect number of input arguments.");
	}
	bool csr = (bool) mxGetScalar(prhs[2]);
	bool bsr = (bool) mxGetScalar(prhs[3]);
	plhs[0] = mxCreateDoubleScalar(core_ptr->is_all_sparse(csr, bsr));
}

template <typename SCALAR, FDevice DEV>
void faust_is_all_dense(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs != 1)
	{
		mexErrMsgTxt("is_all_dense: incorrect number of output arguments.");
	}
	plhs[0] = mxCreateDoubleScalar(core_ptr->is_all_dense());
}
