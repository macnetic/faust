#include "faust_GivensFGFTComplex.h"

#include <list>

#ifndef __GIVENS_FGFT_PARALLEL_COMPLEX__
#define __GIVENS_FGFT_PARALLEL_COMPLEX__

namespace Faust {


	template<typename FPP, Device DEVICE, typename FPP2 = float>
		class GivensFGFTParallelComplex : public GivensFGFTComplex<FPP,DEVICE,FPP2>
	{
		/**
		 * \class Faust::GivensFGFTParallelComplex
		 *
		 * \brief This class implements the parallel version of Givens FGFT algorithm.
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
		 */

		/** Maximum number of rotations per factor
		* the effective number of rotations is the min(t, number_of_pivot_candidates)
		*/
		unsigned int t;
		// Number of rotations already made for the current factor facts[ite]
		unsigned int fact_nrots;
		// current fact nonzero indices (descendantly sorted according to absolute values in L)
		list<pair<int,int>> fact_nz_inds;

		void max_L();

		/**
		 * After selecting a pivot in iteration previous steps, it's necessary to remove them.
		 * That's what this function is responsible for.
		 * The pivots are nonzero elements of L and their indices are sorted in fact_nz_inds (ascendently according to the pivot values).
		 */
		void update_fact_nz_inds();
		/**
		 * This function is responsible to compute all the coefficients (excepting the identity part)
		 * of the current iteration factor (the coefficients are those from the rotation matrices).
		 */
		void loop_update_fact();
		void choose_pivot();
		/**
		 * Computes the coefficients of the last selected rotation matrix to be put later in current iteration factor.
		 */
		void update_fact();
		void update_L(Faust::MatSparse<FPP,Cpu> &L);
		void update_L(Faust::MatDense<FPP,Cpu> &L);
		/**
		 * Constructs the current factor after computing all the coefficients (of rotation matrices) in temporary buffers (see update_fact()).
		 */
		void finish_fact();

		/**
		 * Function pointer to any step of the algorithm.
		 */
		typedef void (GivensFGFTParallelComplex<FPP,DEVICE,FPP2>::*substep_fun)();


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
		GivensFGFTParallelComplex(Faust::MatDense<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity = 0, const double stoppingCritIsError = 0.0,  const bool errIsRel = true);
		GivensFGFTParallelComplex(Faust::MatSparse<FPP,DEVICE>& Lap, int J, int t, unsigned int verbosity = 0, const double stoppingCritIsError = 0.0, const bool errIsRel = true);

	};

#include "faust_GivensFGFTParallelComplex.hpp"

}
#endif
