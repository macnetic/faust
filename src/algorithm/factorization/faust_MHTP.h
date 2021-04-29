#ifndef __MHTP__
#define __MHTP__
namespace Faust
{
	template<typename FPP>
		struct MHTPParams
		{
			bool used;
			StoppingCriterion<FPP> sc;
			bool constant_step_size;
			FPP step_size;
			int palm4msa_period;
			bool updating_lambda;

			MHTPParams() : used(false) {}; //TODO: move in faust_MHTP.hpp
		};
};
#endif
