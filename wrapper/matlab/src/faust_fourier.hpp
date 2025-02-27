#include "class_handle.hpp"
#include "faust_TransformHelper.h"
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#endif
template <typename SCALAR, FDevice DEV>
void faust_fourier(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{

	if(nlhs!=1)
		mexErrMsgTxt("mex dft(): too many left hand side variables.");
	if(nrhs < 4)
		mexErrMsgTxt("mex dft(): wrong number of arguments (must be 4)");
	unsigned int n  = (unsigned int) mxGetScalar(prhs[1]);
	bool norma = (bool) mxGetScalar(prhs[2]);
	bool diag_prod = (bool) mxGetScalar(prhs[3]);
	Faust::TransformHelper<std::complex<double>,DEV> * F = nullptr;
	if(diag_prod)
		F = Faust::TransformHelper<std::complex<double>,DEV>::fourierFaustOpt(n, norma);
	else
		F = Faust::TransformHelper<std::complex<double>,DEV>::fourierFaust(n, norma);

	if(F) //not NULL
	{
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<std::complex<double>,DEV> >(F);
		return;
	}

	// (F == nullptr)
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	double* ptr_out = (double*) mxGetData(plhs[0]);
	ptr_out[0]=(double) 0;
}

