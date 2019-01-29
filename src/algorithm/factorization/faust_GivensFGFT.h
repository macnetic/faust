
#ifndef __GIVENS_FGFT__
#define __GIVENS_FGFT__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include <vector>

namespace Faust {

	template<typename FPP, Device DEVICE>
		class GivensFGFT {

			vector<Faust::MatSparse<FPP,DEVICE>> facts;
			Faust::MatSparse<FPP,DEVICE> D;
			Faust::MatSparse<FPP,DEVICE> C;
			vector<float> err;
			vector<pair<int,int>> coord_choices;
			Faust::MatDense<FPP, DEVICE> L;

			public:
			GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, unsigned int j);
			void next_step();
			void compute_facts();

		};

}
#include "faust_GivensFGFT.hpp"
#endif
