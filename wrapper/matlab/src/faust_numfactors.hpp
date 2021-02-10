#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_numfactors(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs != 1 || nrhs != 2)
	{
		mexErrMsgTxt("numfactors: incorrect number of arguments.");
	}

	faust_unsigned_int nb_fact = core_ptr->size();
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	double* ptr_out = (double*) mxGetData(plhs[0]);
	ptr_out[0]=(double) nb_fact;
}

