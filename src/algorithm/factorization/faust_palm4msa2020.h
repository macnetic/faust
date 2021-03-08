#ifndef __FAUST_PALM4MSA__
#define __FAUST_PALM4MSA__
#include "faust_MatGeneric.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_TransformHelper.h"
#ifdef USE_GPU_MOD
#include "faust_MatGeneric_gpu.h"
#include "faust_MatDense_gpu.h"
#include "faust_MatSparse_gpu.h"
#include "faust_TransformHelper_gpu.h"
#endif
#include "faust_constant.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintMat.h"
#include "faust_StoppingCriterion.h"
#include "faust_prod_opt.h"
#include <functional>

namespace Faust
{

	template <typename FPP, FDevice DEVICE>
		void palm4msa(
				/** input matrix */
				const MatDense<FPP,DEVICE>& A,
				/* input constraints */
				std::vector<ConstraintGeneric*> & constraints,
				/** output factors (if not initialized to size nfacts, set to nfacts identities (with constraints dims) */
				TransformHelper<FPP,DEVICE>& S,
				/** lambda output, intialized from outside */
				FPP& lambda,
				//const unsigned int nites,
				const StoppingCriterion<Real<FPP>>& sc,
				const bool is_update_way_R2L=false,
                const bool use_csr=true,
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
				const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION,
				const bool on_gpu=false,
				const bool is_verbose=false); //TODO: this def of this decl must be updated in faust_palm4msa2020.hpp


	template <typename FPP, FDevice DEVICE>
		void palm4msa2(
				/** input matrix */
				const MatDense<FPP,DEVICE>& A,
				/* input constraints */
				std::vector<ConstraintGeneric*> & constraints,
				/** output factors (if not initialized to size nfacts, set to nfacts identities (with constraints dims) */
				TransformHelper<FPP,DEVICE>& S,
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
				const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION,
				const bool on_gpu=false,
				const bool is_verbose=false,
				/* id argument is useful to identify the palm4msa call.
				 * For example, hierarchical use it to indicate the number of the current factorization.
				 * */
				const int id=0);

	/**
	 * \brief Fill S with nfacts eye() matrices.
	 *
	 * S must be an empty Faust. If S is already composed of factors an exception is raised.
	 *
	 * \param sparse if true eyes are MatSparse, they are MatDense otherwise.
	 *
	 */
	template <typename FPP, FDevice DEVICE>
		void fill_of_eyes(TransformHelper<FPP,DEVICE>& S,
				const unsigned int nfacts,
				const bool sparse,
				const std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims,
				const bool on_gpu=false);
	//TODO: maybe move this function in a utils module (for now it serves only for PALM4MSA)

	// warning: before calling compute_n_apply_grad*() out must be initialized to S[f_id] : the factor to update
	// TODO: ideally compute_n_apply_grad1 has no reason to be kept, compute_n_apply_grad2 is faster (but just in case I decided to keep it for a moment)
	template <typename FPP, FDevice DEVICE>
		void compute_n_apply_grad1(const int f_id, const MatDense<FPP,DEVICE> &A, TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const FPP& lambda, const Real<FPP>& c, MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod, const bool packing_RL);

	template <typename FPP, FDevice DEVICE>
		void compute_n_apply_grad2(const int f_id, const MatDense<FPP,DEVICE> &A, TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const FPP& lambda, const Real<FPP> &c, MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod, const bool packing_RL);

}
#include "faust_palm4msa2020.hpp"
#include "faust_palm4msa2020_2.hpp"
#endif
