#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_norm_2(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs > 5 || nrhs < 3)
		mexErrMsgTxt("norm: Unexpected arguments");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("norm : empty faust core");

	double precision =  0.001;
	int nbr_iter_max=100;
	double norm_faust;
	int flag;

	if(nrhs > 2)
		// the treshold come before
		precision = mxGetScalar(prhs[2]); //mxGetUint32s not available before mat2019
	if(nrhs > 3)
	{
		nbr_iter_max = (int) mxGetScalar(prhs[3]);
		if(nbr_iter_max < 0 || nbr_iter_max > 1000000000)
		{
			mexWarnMsgTxt("Invalid max_num_its value, go back to default value.");
			nbr_iter_max=100;
		}
	}

//	std::cout << "nbr_iter_max:" << nbr_iter_max << " threshold:" << precision << std::endl;
	norm_faust = (double) core_ptr->spectralNorm(nbr_iter_max,precision,flag);

	plhs[0]=mxCreateDoubleScalar(norm_faust);

}

template<typename SCALAR, FDevice DEV>
void get_full_array_and_batch_size(bool& full_array, int& batch_size, const mxArray **prhs, const int nrhs)
{
	// this function doesn't need templates but it's more convenient for compiling/linking
	if(nrhs > 2)
	{
		full_array = mxGetScalar(prhs[2]);
		if(nrhs > 3)
			batch_size = mxGetScalar(prhs[3]);
	}
}

template <typename SCALAR, FDevice DEV>
void faust_norm_1(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	bool full_array = false;
	int batch_size = 1;
	Real<SCALAR> faust_norm;
	get_full_array_and_batch_size<SCALAR, DEV>(full_array, batch_size, prhs, nrhs);
	faust_norm = (double) core_ptr->normL1(full_array, batch_size);
	if(nlhs != 1) mexErrMsgTxt("norm: wrong number of output/left-hand side arguments (it should be one).");
	plhs[0]=mxCreateDoubleScalar(faust_norm);
}

template <typename SCALAR, FDevice DEV>
void faust_norm_inf(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	bool full_array = false;
	int batch_size = 1;
	Real<SCALAR> faust_norm;
	get_full_array_and_batch_size<SCALAR, DEV>(full_array, batch_size, prhs, nrhs);
	faust_norm = (double) core_ptr->normInf(full_array, batch_size);
	if(nlhs != 1) mexErrMsgTxt("norm: wrong number of output/left-hand side arguments (it should be one).");
	plhs[0]=mxCreateDoubleScalar(faust_norm);
}

template <typename SCALAR, FDevice DEV>
void faust_norm_fro(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	bool full_array = false;
	int batch_size = 1;
	Real<SCALAR> faust_norm;
	get_full_array_and_batch_size<SCALAR, DEV>(full_array, batch_size, prhs, nrhs);
	faust_norm = (double) core_ptr->normFro(full_array, batch_size);
	if(nlhs != 1) mexErrMsgTxt("norm: wrong number of output/left-hand side arguments (it should be one).");
	plhs[0]=mxCreateDoubleScalar(faust_norm);
}
