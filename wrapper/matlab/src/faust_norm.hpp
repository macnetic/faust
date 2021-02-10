#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_norm(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs > 5 || nrhs < 3)
		mexErrMsgTxt("norm: Unexpected arguments");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("norm : empty faust core");

	double precision =  0.001;
	double norm_faust;
	int nbr_iter_max=100;
	int flag;
	int ord = (int) mxGetScalar(prhs[2]);

	if(nrhs > 3)
		// the treshold come before
		precision = mxGetScalar(prhs[3]); //mxGetUint32s not available before mat2019
	if(nrhs > 4)
	{
		nbr_iter_max = (int) mxGetScalar(prhs[4]);
		if(nbr_iter_max < 0 || nbr_iter_max > 1000000000)
		{
			mexWarnMsgTxt("Invalid max_num_its value, go back to default value.");
			nbr_iter_max=100;
		}
	}

	if(ord==2)
		norm_faust = (double) core_ptr->spectralNorm(nbr_iter_max,precision,flag);
	else if(ord==1)
		norm_faust = (double) core_ptr->normL1();
	else if(ord == -(1<<(32-1))) //that's the value for 'inf' in matlab
		// a constant from c-c++ would be better
		norm_faust = (double) core_ptr->normInf();
	else
		mexErrMsgTxt("norm : invalid norm order.");

	plhs[0]=mxCreateDoubleScalar(norm_faust);

}

template <typename SCALAR, FDevice DEV>
void faust_norm_fro(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs != 1 || nrhs != 2)
		mexErrMsgTxt("normfro: Unexpected arguments");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("normfro : empty faust core");


	double norm_faust;

	norm_faust = core_ptr->normFro();

	plhs[0]=mxCreateDoubleScalar(norm_faust);

}
