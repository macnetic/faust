#ifndef __FAUST_PALM4MSA__
#define __FAUST_PALM4MSA__
#include "faust_MatGeneric.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_TransformHelper.h"
#include "faust_constant.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintMat.h"
#include "faust_StoppingCriterion.h"
#include <functional>

namespace Faust
{
	template <typename FPP, FDevice DEVICE>
		void palm4msa(
				/** input matrix */
				const Faust::MatDense<FPP,DEVICE>& A,
				/* input constraints */
				std::vector<Faust::ConstraintGeneric*> & constraints,
				/** output factors (if not initialized to size nfacts, set to nfacts identities (with constraints dims) */
				Faust::TransformHelper<FPP,DEVICE>& S,
				/** lambda output, intialized from outside */
				FPP& lambda,
				//const unsigned int nites,
				const StoppingCriterion<Real<FPP>>& sc,
				const bool is_update_way_R2L=false,
                const bool use_csr=true,
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
				const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION);

	template <typename FPP, FDevice DEVICE>
		void palm4msa2(
				/** input matrix */
				const Faust::MatDense<FPP,DEVICE>& A,
				/* input constraints */
				std::vector<Faust::ConstraintGeneric*> & constraints,
				/** output factors (if not initialized to size nfacts, set to nfacts identities (with constraints dims) */
				Faust::TransformHelper<FPP,DEVICE>& S,
				/** lambda output, intialized from outside */
				Real<FPP>& lambda,
				//const unsigned int nites,
				const StoppingCriterion<Real<FPP>>& sc,
				const bool is_update_way_R2L=false,
                const bool use_csr=true,
				const bool packing_RL=true,
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
				const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION);

}
#include "faust_palm4msa2020.hpp"
#endif
