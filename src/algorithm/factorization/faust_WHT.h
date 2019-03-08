#ifndef FAUST_WHT
#define FAUST_WHT
#include "faust_MatGeneric.h"
#include <vector>

namespace Faust {
	/**
	 * \brief Fast Walsh-Hadamard Transform.
	 */
	template<typename FPP>
	void wht_factors(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  factors);

}
#include "faust_WHT.hpp"
#endif
