#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_hadamard(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	if(nlhs!=1)
		mexErrMsgTxt("mex wht(): too many left hand side variables.");
	if(nrhs < 3)
		mexErrMsgTxt("mex wht(): wrong number of arguments (must be 3)");
	unsigned int n  = (unsigned int) mxGetScalar(prhs[1]);
	bool norma = (bool) mxGetScalar(prhs[2]);

	Faust::TransformHelper<SCALAR,DEV>* H = Faust::TransformHelper<SCALAR,DEV>::hadamardFaust(n, norma);
	//H->display();
	if(H) //not NULL
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(H);
	else {
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}

}

