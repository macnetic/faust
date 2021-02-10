#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_size(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV>>(prhs[1]);
	const size_t SIZE_B1 = core_ptr->getNbRow();
	const size_t SIZE_B2 = core_ptr->getNbCol();
	const mwSize dims[2]={1,2};
	plhs[0]=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
	double* ptr_out = (double*) mxGetData(plhs[0]);
	ptr_out[0]=(double) SIZE_B1;
	ptr_out[1]=(double) SIZE_B2;

	return;
}
