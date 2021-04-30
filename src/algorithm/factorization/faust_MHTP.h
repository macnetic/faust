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
};
#include "faust_MHTP.hpp"
#endif
