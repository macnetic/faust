#ifndef __FAUST_HIERARCHICAL__
#define __FAUST_HIERARCHICAL__
#include "faust_palm4msa2020.h"
namespace Faust
{
	template<typename FPP, Device DEVICE>
	Faust::TransformHelper<FPP,DEVICE>* hierarchical(const Faust::MatDense<FPP,DEVICE>& A,
				const int nites,
				std::vector<const Faust::ConstraintGeneric*> & fac_constraints,
				std::vector<const Faust::ConstraintGeneric*> & res_constraints,
				Real<FPP>& lambda,
				const bool is_update_way_R2L=false, const bool is_fact_side_left=false,
				const bool use_csr=true, const bool packing_RL=true,
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER);
}
#include "faust_hierarchical.hpp"
#endif
