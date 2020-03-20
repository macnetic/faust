#include "faust_GivensFGFTComplex.h"
#include "faust_GivensFGFTParallelGen.h"

#include <list>

#ifndef __GIVENS_FGFT_PARALLEL_COMPLEX__
#define __GIVENS_FGFT_PARALLEL_COMPLEX__

namespace Faust {


	template<typename FPP, Device DEVICE, typename FPP2 = float>
		class GivensFGFTParallelComplex : public GivensFGFTComplex<FPP,DEVICE,FPP2>, public GivensFGFTParallelGen<typename FPP::value_type, DEVICE, FPP2, FPP>
	{
		/**
		 * \class Faust::GivensFGFTParallelComplex
		 *
		 * \brief This class implements the parallel version of Givens FGFT algorithm (for the complex case).
		 *
		 * This variant of the parent class algorithm consists mainly to put t 2D rotation matrices in each iteration factor S (i.e. facts[ite]) when the basis version puts only a single rotation matrix into L.
		 *
		 * This algorithm is based on the classical Jacobi eigenvalues algorithm.
		 *
		 *  References:
		 *
		 *  [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
		 *    graph Fourier transforms via multi-layer sparse approximations",
		 *    submitted to IEEE Transactions on Signal and Information Processing
		 *    over Networks.
		 *    <https://hal.inria.fr/hal-01416110>
		 *
		 *  Complementary reference for the complex case: https://en.wikipedia.org/wiki/Jacobi_method_for_complex_Hermitian_matrices
		 *
		 */


		/**
		 * Computes the coefficients of the last selected rotation matrix to be put later in current iteration factor.
		 */
		void update_fact();
		void update_L(Faust::MatDense<FPP,Cpu> &L);

		/**
		 * Function pointer to any step of the algorithm.
		 */
		typedef void (GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::*substep_fun)();


		void init_fact_nz_inds_sort_func();

		public:

		void next_step();

		/**
		 * Constructor.
		 *
		 *
		 * \param Lap The Laplacian matrix to approximate/diagonalize.
		 * \param J round(J/t) is the number of factors to approximate the Fourier matrix (as a product of factors).
		 * \param t the maximum number of 2D rotation matrices (Givens matrices) to insert in each factor. The effective number can be less than t if there is not enough pivot candidates left in the matrix L of the current iteration.
		 * \param stoppingError defines a stopping criterion based on error (absolute relative error).
		 *
		 */
		GivensFGFTParallelComplex(Faust::MatDense<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity = 0, const double stoppingCritIsError = 0.0,  const bool errIsRel = true, const bool enable_large_Faust = false);
		GivensFGFTParallelComplex(Faust::MatSparse<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity = 0, const double stoppingCritIsError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false);
	};

#include "faust_GivensFGFTParallelComplex.hpp"

}
#endif
