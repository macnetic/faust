#ifndef __FAUST_HIERARCHICAL__
#define __FAUST_HIERARCHICAL__
#include "faust_palm4msa2020.h"
#include "faust_Params.h"
namespace Faust
{
	template<typename FPP, FDevice DEVICE>
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

	template<typename FPP, FDevice DEVICE>
		Faust::TransformHelper<FPP,DEVICE>* hierarchical(const Faust::MatDense<FPP,DEVICE>&  A,
				Faust::Params<FPP,DEVICE, Real<FPP>> &p,
				Real<FPP>& lambda, const bool compute_2norm_on_array);
}

#include "faust_hierarchical.hpp"
#endif
