#include "class_handle.hpp"
#include "faust_TransformHelper.h"
#include "faust2Mx.h"
template <typename SCALAR, FDevice DEV>
void faust_full(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	if (nlhs != 1 || nrhs != 2)
		mexErrMsgTxt("full: Unexpected arguments");
	if(core_ptr->size() == 0)
		mexErrMsgTxt("full : empty faust core");

	faust_unsigned_int nbRowOp,nbColOp;
	const size_t SIZE_B1 = core_ptr->getNbRow();
	const size_t SIZE_B2 = core_ptr->getNbCol();
	Faust::MatDense<SCALAR,Cpu> prod;
	core_ptr->get_product(prod);

	const mwSize dims[2]={SIZE_B1,SIZE_B2};
	plhs[0] = FaustMat2mxArray(prod);
}

