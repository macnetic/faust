#ifndef __MHTP__
#define __MHTP__
#include <string>
#include "faust_StoppingCriterion.h"
#include "faust_ConstraintGeneric.h"
namespace Faust
{
	/**
	 * \brief This class represents the set of parameters used for the MHTP (Multilinear Hard Thresholding Pursuit) algorithm used optionally in PALM4MSA (2020 implementation only).
	 *
	 * \param used: a boolean to indicate if the MHTP algorithm will be launched (and this set of parameters used). It makes sense in PALM4MSA which can launch MHTP or no.
	 * \param sc: The stopping criterion used to determine the number of iterations of the MHTP algorithm/pass (in PALM4MSA).
	 * \param constant_step_size: true to set a constant step size in the gradient descent of MHTP. Otherwise (if constant_step_size is false) then the gradient is recomputed at each iteration according to the lipschitz coefficient.
	 * \param step_size: the step size used if constant_step_size is true.
	 * \param updating_lambda: true to recompute the scaling factor lambda of the Faust at each iteration of MHTP (if false lamdba is updated only at each iteration of PALM4MSA).
	 */
	template<typename FPP>
		struct MHTPParams
		{
			bool used;
			StoppingCriterion<FPP> sc;
			bool constant_step_size;
			FPP step_size;
			int palm4msa_period;
			bool updating_lambda;
			MHTPParams();
			std::string to_string() const;
		};

	/**
	 * \brief This function performs a MHTP pass on one factor of S for as many iterations as defined in mhtp_params.
	 *
	 * This function is intended to be called from a PALM4MSA implementation (especially PALM4MSA 2020 implementation).
	 *
	 * \param mhtp_params: the MHTPParams to configure the MHTP algorithm.
	 * \param A: The matrix being approximated.
	 * \param A_H: The transconjugate of the matrix being approximated (it avoids to recompute for each call).
	 * \param S: The Faust (multilayer matrix) approximate of A to refine.
	 * \param f_id: The factor on which this MHTP pass will work.
	 * \param pL: The left-hand side Faust of S[f_id].
	 * \param pR: The right-hand side Faust of S[f_id].
	 * \param is_verbose: if true the function will be verbose, otherwise no information is displayed.
	 * \param constraint: the constraint to apply on S[f_id]
	 * \param norm2_max_iter: the maximum number of iterations (of power iteration algorithm) for the 2-norm computed.
	 * \param norm2_threshold: the threshold/precision according to compute the 2-norm.
	 * \param norm2_duration: measurement of the time passed computing 2-norms in this function (if norm2_duration is not zero at call time, the time is added instead of overridden).
	 * \param fgrad_duration: measurment of the time passed computing/applying the gradient in this function.
	 * \param sc: The StoppingCriterion used in the gradient descent used to update the factor (as in PALM4MSA).
	 * \param error: after this function execution error will be equal to norm(A-S) (frobenius norm).
	 * \param use_csr: true if S is only sparse factors, false if they are all dense.
	 * \param prod_mod: the method used when computing the full matrix of S (see TransformHelper::get_product() and FaustMulMode).
	 * \param c: define step size of the gradient descent (which is set by the function to 1/mhtp_params.step_size if mhtp_params.constant_step_size == true, else it is computed according to the lipschitz coefficient).
	 * \param lambda: the scaling factor of S (lambda*S approx A) which is modified internally if mhtp_params.updating_lambda is true.
	 *
	 * Reference: https://hal.inria.fr/hal-03132013/document
	 */
	template<typename FPP, FDevice DEVICE>
		void perform_MHTP(
				const MHTPParams<FPP>& mhtp_params,
				const Faust::MatDense<FPP,DEVICE>& A,
				const Faust::MatDense<FPP,DEVICE>& A_H,
				Faust::TransformHelper<FPP,DEVICE>& S,
				int f_id,
				std::vector<TransformHelper<FPP,DEVICE>*> &pL,
				std::vector<TransformHelper<FPP,DEVICE>*> &pR,
				const bool is_verbose,
				const Faust::ConstraintGeneric &constraint,
				const int norm2_max_iter,
				const Real<FPP>& norm2_threshold,
				std::chrono::duration<double>& norm2_duration,
				std::chrono::duration<double>& fgrad_duration,
				const StoppingCriterion<Real<FPP>>& sc,
				Real<FPP> &error,
				const bool use_csr,
				const int prod_mod,
				Real<FPP> &c,
				Real<FPP>& lambda);

};
#include "faust_MHTP.hpp"
#endif
