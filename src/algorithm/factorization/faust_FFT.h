#ifndef FAUST_FFT
#define FAUST_FFT
#include "faust_MatGeneric.h"
#include <vector>

namespace Faust {
	/**
	 * \brief Cooley-tukey FFT algorithm.
	 * Ref.: http://www.cs.cornell.edu/~bindel/class/cs5220-s10/slides/FFT.pdf
	 */
	template<typename FPP>
	void fft_factors(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  v);

}
#include "faust_FFT.hpp"
#endif
