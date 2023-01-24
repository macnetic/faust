#include "faust2Mx.h"
#include "class_handle.hpp"
#include "faust_TransformHelper.h"

/**
 * These functions allow to convert a MatBSR factor of a Faust::TransformHelper to MatSparse, as the bsr matrices are not supported in Matlab.
 */

template <typename SCALAR, FDevice DEV>
void faust_factors(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs > 1 || nrhs != 3)
	{
		mexErrMsgTxt("factors : incorrect number of arguments.");
	}
	int id = (int) (mxGetScalar(prhs[2])-1);
	auto type = core_ptr->get_fact_type(id);

	if(core_ptr->is_fact_sparse(id))
		plhs[0] = transformFact2SparseMxArray(id,core_ptr);
//
//		plhs[0] = FaustSpMat2mxArray(*dynamic_cast<const Faust::MatSparse<SCALAR,Cpu>*>(core_ptr->get_gen_fact(id)));
	else if(core_ptr->is_fact_dense(id))
		plhs[0] = transformFact2FullMxArray(id,core_ptr);
//		plhs[0] = FaustMat2mxArray(*dynamic_cast<const Faust::MatDense<SCALAR,Cpu>*>(core_ptr->get_gen_fact(id)));
	else if(core_ptr->is_fact_bsr(id))
	{
		plhs[0] = bsr_mat_to_sp_mat(id, core_ptr);
	}
	else if(type == 4)
		// MatButterfly
		plhs[0] = butterfly_mat_to_sp_mat(id, core_ptr);
	else if(type == 5)
		// MatPerm
		plhs[0] = perm_mat_to_sp_mat(id, core_ptr);
	else
		mexErrMsgTxt("Unhandled type of matrix");

}

template <typename SCALAR>
mxArray* bsr_mat_to_sp_mat(int id, Faust::TransformHelper<SCALAR,Cpu>* core_ptr)
{
	auto bsr_mat = dynamic_cast<const MatBSR<SCALAR, Cpu>*>(core_ptr->get_gen_fact(id));
	auto sp_mat = bsr_mat->to_sparse();
	return FaustSpMat2mxArray(sp_mat);
}


template <typename SCALAR>
mxArray* butterfly_mat_to_sp_mat(int id, Faust::TransformHelper<SCALAR,Cpu>* core_ptr)
{
	auto bf_mat = dynamic_cast<const MatButterfly<SCALAR, Cpu>*>(core_ptr->get_gen_fact(id));
	auto sp_mat = bf_mat->toMatSparse();
	return FaustSpMat2mxArray(sp_mat);
}


template <typename SCALAR>
mxArray* perm_mat_to_sp_mat(int id, Faust::TransformHelper<SCALAR,Cpu>* core_ptr)
{
	auto p_mat = dynamic_cast<const MatPerm<SCALAR, Cpu>*>(core_ptr->get_gen_fact(id));
	auto sp_mat = p_mat->toMatSparse();
	return FaustSpMat2mxArray(sp_mat);
}

#ifdef USE_GPU_MOD
template <typename SCALAR>
mxArray* bsr_mat_to_sp_mat(int id, Faust::TransformHelper<SCALAR, GPU2>* core_ptr)
{
	//TODO: nothing to do now because MatBSR aren't supported on GPU2, implement it later
	return nullptr;
}
#endif
