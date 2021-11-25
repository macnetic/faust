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
#include "faust_linear_algebra_gpu.h"
#include "faust_TransformHelper_gpu.h"
#endif
#include "faust_Params.h"
#include "faust_constant.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintMat.h"
#include "faust_StoppingCriterion.h"
#include "faust_prod_opt.h"
#include "faust_MHTP.h"
#include <functional>
#include <cstdlib>
#include <cmath>

#define PALM4MSA2020_VERBOSE_CALC_ERR_ITE_PERIOD 1 // the period according to the relative error is computed and displayed in palm4msa2
// this constant is overriden if the variable environment VERBOSE_CALC_ERR_ITE_PERIOD exists
#define set_calc_err_ite_period() \
			auto ite_period = PALM4MSA2020_VERBOSE_CALC_ERR_ITE_PERIOD; \
			auto env_var = getenv("PALM4MSA2020_VERBOSE_CALC_ERR_ITE_PERIOD"); \
			if (nullptr != env_var) \
				ite_period = atoi(env_var);

#define LIPSCHITZ_MULTIPLICATOR 1.001

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
				const FactorsFormat factors_format=AllDynamic,
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
				const bool constant_step_size=false, const Real<FPP> step_size=FAUST_PRECISION,
				const bool on_gpu=false,
				const bool is_verbose=false); //TODO: the def of this decl must be updated in faust_palm4msa2020.hpp


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
				const FactorsFormat factors_format=AllDynamic,
				const bool packing_RL=true,
				const bool no_normalization=false,
				const bool no_lambda=false,
				const MHTPParams<Real<FPP>> mhtp_params=MHTPParams<Real<FPP>>(),
				const bool compute_2norm_on_array=false,
				const Real<FPP> norm2_threshold=FAUST_PRECISION,
				const unsigned int norm2_max_iter=FAUST_NORM2_MAX_ITER,
				bool constant_step_size=false, Real<FPP> step_size=FAUST_PRECISION,
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
		void compute_n_apply_grad1(const int f_id, const MatDense<FPP,DEVICE> &A, TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const bool packing_RL, const Real<FPP>& lambda, const Real<FPP>& c, MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod);

	template <typename FPP, FDevice DEVICE>
		void compute_n_apply_grad2(const int f_id, const MatDense<FPP,DEVICE> &A, TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR,  const bool packing_RL, const Real<FPP>& lambda, const Real<FPP> &c, MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod);

	template<typename FPP, FDevice DEVICE>
		Real<FPP> calc_rel_err(const TransformHelper<FPP,DEVICE>& S, const MatDense<FPP,DEVICE> &A, const Real<FPP> &lambda=1, const Real<FPP>* A_norm=nullptr);

	/**
	 * \brief This function performs the (scaling factor) lambda update of the PALM4MSA algorithm (palm4msa2).
	 *
	 * \param S: the Faust being refined by the PALM4MSA algorithm (palm4msa2).
	 * \param A_H: the transconjugate of the matrix A for which S is an approximate.
	 * \param lambda: the output of the lambda computed by the function.
	 */
	template<typename FPP, FDevice DEVICE>
		void update_lambda(Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const MatDense<FPP, DEVICE> &A_H, Real<FPP>& lambda, bool no_lambda_error=false);

	template<typename FPP, FDevice DEVICE>
		void update_fact(
				Faust::MatGeneric<FPP,DEVICE>* cur_fac,
				int f_id,
				const Faust::MatDense<FPP,DEVICE>& A,
				Faust::TransformHelper<FPP,DEVICE>& S,
				std::vector<TransformHelper<FPP,DEVICE>*> &pL,
				std::vector<TransformHelper<FPP,DEVICE>*> &pR,
				const bool packing_RL,
				const bool is_verbose,
				const Faust::ConstraintGeneric &constraints,
				const int norm2_max_iter,
				const Real<FPP>& norm2_threshold,
				std::chrono::duration<double>& spectral_duration,
				std::chrono::duration<double>& fgrad_duration,
				const bool constant_step_size,
				const Real<FPP> step_size,
				const StoppingCriterion<Real<FPP>>& sc,
				Real<FPP> &error,
				const FactorsFormat factors_format,
				const int prod_mod,
				Real<FPP> &c,
				const Real<FPP>& lambda,
				bool use_grad1=false);

	/**
	 * These two functions compute the mat 2-norm using double precision. Intended to be used only when FPP is float.
	 * It is a workaround to the 2-norm float computation (NaN result) which occur on certain Fausts (cf. gitlab issue #236).
	 */
	template<typename FPP>
		Real<FPP> compute_double_spectralNorm(MatDense<FPP, Cpu>& mat, int norm2_max_iter, double norm2_threshold);
	template<typename FPP>
		Real<FPP> compute_double_spectralNorm(MatDense<FPP, GPU2>& mat, int norm2_max_iter, double norm2_threshold);


}
#include "faust_palm4msa2020.hpp"
#include "faust_palm4msa2020_2.hpp"
#endif
