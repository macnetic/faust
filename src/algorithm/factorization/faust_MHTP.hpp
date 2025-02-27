#include "faust_palm4msa2020.h"
namespace Faust
{
	template<typename FPP>
		MHTPParams<FPP>::MHTPParams() : used(false),
		constant_step_size(false),
		step_size(1e-3),
		palm4msa_period(1000),
		updating_lambda(true),
		sc(StoppingCriterion<FPP>(50)){};

	template<typename FPP>
	std::string MHTPParams<FPP>::to_string() const
	{
		auto sc_str = sc.to_string();
		std::string str = "MHTPParams (START):";
		str += "\r\n";
		str += "StoppingCriterion:";
		str += "\r\n === \r\n";
		str += sc_str;
		str += "\r\n === \r\n";
		str += "constant_step_size: ";
		str += std::to_string(constant_step_size);
		str += "\r\n";
		str += "step_size: ";
		str += std::to_string(step_size);
		str += "\r\n";
		str += "palm4msa_period: ";
		str += std::to_string(palm4msa_period);
		str += "\r\n";
		str += "updating_lambda: ";
		str += std::to_string(updating_lambda);
		str += "\r\n";
		str += "MHTPParams END.";
		return str;
	}

	template<typename FPP, FDevice DEVICE>
        void perform_MHTP(
                const int it,
                const MHTPParams<Real<FPP>>& mhtp_params,
                const Faust::MatDense<FPP,DEVICE>& A,
                const Faust::MatDense<FPP,DEVICE>& A_H,
                Faust::TransformHelper<FPP,DEVICE>& S,
                const int f_id,
                const bool is_update_way_R2L,
                std::vector<TransformHelper<FPP,DEVICE>*> &pL,
                std::vector<TransformHelper<FPP,DEVICE>*> &pR,
                const bool packing_RL,
                const bool is_verbose,
                const Faust::ConstraintGeneric &constraint,
                const int norm2_max_iter,
                const Real<FPP>& norm2_threshold,
                std::chrono::duration<double>& norm2_duration,
                std::chrono::duration<double>& fgrad_duration,
                const StoppingCriterion<Real<FPP>>& sc,
                Real<FPP> &error,
                const FactorsFormat factors_format,
                const int prod_mod,
                Real<FPP> &c,
                Real<FPP>& lambda)
				{

					Faust::MatGeneric<FPP,DEVICE>* cur_fac;
					// set the factor to zero
					S.get_gen_fact_nonconst(f_id)->setZeros();
					if(is_verbose)
						std::cout << "Starting a MHTP pass ("<< mhtp_params.sc.get_crit() <<" iterations) for factor #" << f_id << std::endl;
					int j = 0;
					while(mhtp_params.sc.do_continue(j)) // TODO: what about the error stop criterion?
					{
						cur_fac = S.get_gen_fact_nonconst(f_id);
						update_fact(it, cur_fac, f_id, is_update_way_R2L, A, S, pL, pR, packing_RL,
								is_verbose, constraint,
								norm2_max_iter, norm2_threshold, norm2_duration, fgrad_duration,
								mhtp_params.constant_step_size, mhtp_params.step_size,
								sc, error, factors_format, prod_mod, c, lambda);
						if(mhtp_params.updating_lambda)
							update_lambda(S, pL, pR, A_H, lambda);
						j++;
					}
					if(is_verbose)
						std::cout << "The MHTP pass has ended" << std::endl;
				}
}


