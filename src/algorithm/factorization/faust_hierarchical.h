#ifndef __FAUST_HIERARCHICAL__
#define __FAUST_HIERARCHICAL__
#include "faust_palm4msa2020.h"
namespace Faust
{
	template<typename FPP, Device DEVICE>
	Faust::TransformHelper<FPP,DEVICE>* hierarchical(const Faust::MatDense<FPP,DEVICE>& A,
//				const int nites,
			std::vector<StoppingCriterion<Real<FPP>>>& sc,
			std::vector<const Faust::ConstraintGeneric*> & fac_constraints,
			std::vector<const Faust::ConstraintGeneric*> & res_constraints,
			Real<FPP>& lambda,
			const bool is_update_way_R2L=false, const bool is_fact_side_left=false,
			const bool use_csr=true, const bool packing_RL=true,
			const bool compute_2norm_on_array=false,
			const Real<FPP> norm2_threshold=FAUST_PRECISION,
			const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
			const bool is_verbose=false,
			const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION);
}

// this macro is only for hierarchical
#define DISPLAY_PARAMS() \
if(is_verbose)\
{\
	std::cout << "Faust::hierarchical2020 parameters:" << endl;\
	std::cout << "nfacts: " << fac_constraints.size()+1 << endl;\
	std::cout << "verbose: " << is_verbose << endl;\
	std::cout << "updateway R2L: " << is_update_way_R2L << endl;\
	std::cout << "init_lambda: " << lambda << endl;\
	std::cout << "is_fact_side_left: " << is_fact_side_left << endl;\
	std::cout << "stop_crit_2facts: " << sc[0].get_crit() << endl;\
	std::cout << "stop_crit_global: " << sc[1].get_crit() << endl;\
	std::cout << "norm2_threshold: " << norm2_threshold << endl;\
	std::cout << "norm2_max_iter: " << norm2_max_iter << endl;\
	std::cout << "packing_RL:" << packing_RL << endl;\
	std::cout << "use_csr:" << use_csr << endl;\
	std::cout << "constant_step_size:" << constant_step_size << endl;\
	std::cout << "step_size:" << step_size << endl;\
	std::cout << "constraints: " << fac_constraints.size() << endl;\
	cout << "FACTORS:" << endl;\
	for(int i=0;i<fac_constraints.size();i++)\
	fac_constraints[i]->Display();\
	cout << "RESIDUUMS:" << endl;\
	for(int i=0;i<res_constraints.size();i++)\
	res_constraints[i]->Display();\
}
//		std::cout << "is_constant_step_size: " << is_constant_step_size << endl;\
//		std::cout << "step_size: " << step_size << endl;\


#include "faust_hierarchical.hpp"
#endif
