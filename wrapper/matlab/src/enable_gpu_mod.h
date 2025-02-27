#ifndef ENABLE_GPU_MOD
#define ENABLE_GPU_MOD
// this function is defined in a header because it is not intended to be used in several mexs
void enable_gpu_mod(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	char libpath[1024];
	char backend[32];
	bool silent;
	if(nrhs < 4)
		mexErrMsgTxt("enable_gpu_mod received a bad number of arguments.");

	mxGetString(prhs[1], libpath, sizeof(libpath));
	mxGetString(prhs[2], backend, sizeof(backend));
	silent = (bool) mxGetScalar(prhs[3]);
	// backend is not yet used (only CUDA is available)
	void *ptr = nullptr;
#ifdef USE_GPU_MOD
	ptr = Faust::enable_gpu_mod(libpath, silent);
#endif
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	double* ptr_out = (double*) mxGetData(plhs[0]);
	ptr_out[0] = (double)((unsigned long long)ptr);
}

void _is_gpu_mod_enabled(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	bool enabled;
#ifdef USE_GPU_MOD
	enabled = Faust::is_gpu_mod_enabled();
	plhs[0] = mxCreateDoubleScalar(enabled);
#else
	enabled = false;
	plhs[0] = mxCreateDoubleScalar(enabled);
#endif
}
#endif
