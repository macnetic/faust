#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_power_ite(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs < 4)
		mexErrMsgTxt("power_ite: need 4 input arguments, 1 output argument.");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("norm : empty faust core");

	double precision =  0.001;
	double eigen_faust; // TODO: support of complex Faust (see the gitlab issue)
	int nbr_iter_max=100;
	int flag;

	precision = mxGetScalar(prhs[2]);
	nbr_iter_max = (int) mxGetScalar(prhs[3]);
	if(nbr_iter_max < 0 || nbr_iter_max > 1000000000)
	{
		mexWarnMsgTxt("Invalid max_num_its value, go back to default value.");
		nbr_iter_max=100;
	}

	auto lambda = core_ptr->power_iteration(nbr_iter_max, precision, flag);

	plhs[0] = mxCreateDoubleScalar(Faust::fabs(lambda)); // fabs to compile with complex even if not used yet (see above todo)
}
