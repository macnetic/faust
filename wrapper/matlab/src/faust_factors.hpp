#include "faust2Mx.h"
#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_factors(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs > 1 || nrhs != 3)
	{
		mexErrMsgTxt("factors : incorrect number of arguments.");
	}
	int id = (int) (mxGetScalar(prhs[2])-1);

	if(core_ptr->is_fact_sparse(id))
		plhs[0] = transformFact2SparseMxArray(id,core_ptr);
//
//		plhs[0] = FaustSpMat2mxArray(*dynamic_cast<const Faust::MatSparse<SCALAR,Cpu>*>(core_ptr->get_gen_fact(id)));
	else
		plhs[0] = transformFact2FullMxArray(id,core_ptr);
//		plhs[0] = FaustMat2mxArray(*dynamic_cast<const Faust::MatDense<SCALAR,Cpu>*>(core_ptr->get_gen_fact(id)));

}

