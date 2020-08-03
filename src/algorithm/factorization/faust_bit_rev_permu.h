#ifndef __FAUST_BIT_REV_PERMU__
#define __FAUST_BIT_REV_PERMU__
namespace Faust {
	/**
	 * \brief This function is only utility for faust_FFT.h(pp).
	 *
	 *  Init v to (0,..,2^n-1).
	 *  And proceeds with a bit reversal permutation.
	 *  E.g. : if n = 3, v=(0, 1, 2, 3, 4, 5, 6, 7) and becomes (0, 4, 2, 6, 1, 5, 3, 7).
	 */
	void bit_rev_permu(unsigned int n, unsigned int* v, const bool initv=true);
}
#endif
