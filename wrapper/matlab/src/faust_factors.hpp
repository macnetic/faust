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
		mexErrMsgTxt("factors : incorrect number of arguments.");

	auto ids = static_cast<unsigned long int*>(mxGetData(prhs[2]));
	auto n_ids = mxGetNumberOfElements(prhs[2]);

	if(n_ids == 1)
	{ // asked a single factor
		int id = ids[0];
		auto type = core_ptr->get_fact_type(id);

		if(core_ptr->is_fact_sparse(id))
			plhs[0] = transformFact2SparseMxArray(id, core_ptr);
		else if(core_ptr->is_fact_dense(id))
			plhs[0] = transformFact2FullMxArray(id, core_ptr);
		else if(core_ptr->is_fact_bsr(id))
			plhs[0] = bsr_mat_to_sp_mat(id, core_ptr);
		else if(type == 4)
			// MatButterfly
			plhs[0] = butterfly_mat_to_sp_mat(id, core_ptr);
		else if(type == 5)
			// MatPerm
			plhs[0] = perm_mat_to_sp_mat(id, core_ptr);
		else
			mexErrMsgTxt("Unhandled type of matrix");
	}
	else
	{ // asked a group of factors // return a Faust
		Faust::TransformHelper<SCALAR, DEV>* th = core_ptr->factors(ids, n_ids);
		plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
	}

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
