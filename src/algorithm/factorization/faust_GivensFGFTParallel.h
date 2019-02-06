#include "faust_GivensFGFT.h"

#include <list>

#ifndef __GIVENS_FGFT_PARALLEL__
#define __GIVENS_FGFT_PARALLEL__

namespace Faust {


	template<typename FPP, Device DEVICE, typename FPP2 = float>
		class GivensFGFTParallel : public GivensFGFT<FPP,DEVICE,FPP2>
	{

		/** Maximum number of rotations per factor
		* the effective number of rotations is the min(t, number_of_pivot_candidates)
		*/
		unsigned int t;
		// Number of rotations already made for the current factor facts[ite]
		unsigned int fact_nrots;
		// current fact nonzero indices (descendantly sorted)
		list<pair<int,int>> fact_nz_inds;

		void max_L();
		void update_fact_nz_inds();
		void loop_update_fact();
		void choose_pivot();
		void update_fact();
		void finish_fact();

		/**
		 * Function pointer to any step of the algorithm.
		 */
		typedef void (GivensFGFTParallel<FPP,DEVICE,FPP2>::*substep_fun)();


		public:

		void next_step();

		GivensFGFTParallel(Faust::MatDense<FPP,DEVICE>& Lap, int J, int t);

	};

#include "faust_GivensFGFTParallel.hpp"


}
#endif
