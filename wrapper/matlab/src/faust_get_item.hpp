#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_get_item(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs != 4)
		mexErrMsgTxt("get_item: Unexpected number of arguments");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("get_item: empty faust core");

	unsigned long int i, j;
	i = (unsigned long int) mxGetScalar(prhs[2]);
	j = (unsigned long int) mxGetScalar(prhs[3]);
	SCALAR item;

	item = core_ptr->get_item(i,j);

	Faust::MatDense<SCALAR, Cpu> mat(1, 1, &item);

	plhs[0] = FaustMat2mxArray(mat);
}
