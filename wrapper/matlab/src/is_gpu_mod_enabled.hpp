template<FDevice DEV>
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
