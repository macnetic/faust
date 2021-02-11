#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_fourier(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{

	if(nlhs!=1)
		mexErrMsgTxt("mex dft(): too many left hand side variables.");
	if(nrhs < 3)
		mexErrMsgTxt("mex dft(): wrong number of arguments (must be 3)");
	unsigned int n  = (unsigned int) mxGetScalar(prhs[1]);
	bool norma = (bool) mxGetScalar(prhs[2]);

	Faust::TransformHelper<complex<Real<SCALAR>>,DEV>* F = Faust::TransformHelper<complex<Real<SCALAR>>,DEV>::fourierFaust(n, norma);
	if(F) //not NULL
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<complex<Real<SCALAR>>,DEV> >(F);
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}
}

