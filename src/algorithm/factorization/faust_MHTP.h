#ifndef __MHTP__
#define __MHTP__
#include <string>
namespace Faust
{
	/**
	 * \brief This class represents the set of parameters used for the MHTP (Multilinear Hard Thresholding Pursuit) algorithm used optionally in PALM4MSA (2020 implementation only).
	 *
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
	 * \brief This function performs the Multilinear Hard Tresholding Pursuit) for as number of iterations as defined in mhtp_params.
	 *
	 * Reference: https://hal.inria.fr/hal-03132013/document
	 */
	template<typename FPP, FDevice DEVICE>
		void perform_MHTP(
				const MHTPParams<FPP>& mhtp_params,
				Faust::MatGeneric<FPP,DEVICE>* cur_fac,
				int f_id,
				const Faust::MatDense<FPP,DEVICE>& A,
				const Faust::MatDense<FPP,DEVICE>& A_H,
				Faust::TransformHelper<FPP,DEVICE>& S,
				std::vector<TransformHelper<FPP,DEVICE>*> &pL,
				std::vector<TransformHelper<FPP,DEVICE>*> &pR,
				const bool is_verbose,
				std::vector<Faust::ConstraintGeneric*> & constraints,
				const int norm2_max_iter,
				const Real<FPP>& norm2_threshold,
				std::chrono::duration<double>& spectral_duration,
				std::chrono::duration<double>& fgrad_duration,
				const StoppingCriterion<Real<FPP>>& sc,
				Real<FPP> &error,
				const bool use_csr,
				const bool packing_RL,
				const int prod_mod,
				Real<FPP> &c,
				Real<FPP>& lambda);

};
#include "faust_MHTP.hpp"
#endif
