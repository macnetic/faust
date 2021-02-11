#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_eye(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	if(nlhs != 1)
		mexErrMsgTxt("eye(): too many left hand side variables.");
	unsigned int n, m;
	n = mxGetScalar(prhs[1]);
	m = mxGetScalar(prhs[2]);
	Faust::TransformHelper<SCALAR,DEV>* eye = Faust::TransformHelper<SCALAR,DEV>::eyeFaust(n, m);
	if(eye) // not NULL
		plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV>>(eye);
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}

}

