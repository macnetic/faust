

#ifndef __GIVENS_FGFT_PARALLEL_GEN__
#define __GIVENS_FGFT_PARALLEL_GEN__
#include "faust_GivensFGFT.h"
#include "faust_MatSparse.h"
#include <list>
#include <functional>
namespace Faust {


	template<typename FPP, Device DEVICE, typename FPP2 = float, typename FPP4 = FPP>
		class GivensFGFTParallelGen
		{
			/**
			 * \class Faust::GivensFGFTParallelGen
			 *
			 * \brief This class represents the parallel version of Givens FGFT algorithm (for concrete implementations see subclasses: GivensFGFTParallelComplex and GivensFGFTParallel).
			 *
			 * This variant of the parent class algorithm consists mainly to put t 2D rotation matrices in each iteration factor S (i.e. facts[ite]) when the basis version puts only a single rotation matrix into L.  //		 * //		 * This algorithm is based on the classical Jacobi eigenvalues algorithm.  //		 * //		 *  References:
			 *
			 *  [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
			 *    graph Fourier transforms via multi-layer sparse approximations",
			 *    submitted to IEEE Transactions on Signal and Information Processing
			 *    over Networks.
			 *    <https://hal.inria.fr/hal-01416110>
			 *
			 */

			Faust::GivensFGFTGen<FPP, DEVICE, FPP2, FPP4> & alg;

			protected:
				/** Maximum number of rotations per factor
				 * the effective number of rotations is the min(t, number_of_pivot_candidates)
				 */
				unsigned int t;
				// Number of rotations already made for the current factor facts[ite]
				unsigned int fact_nrots;
				// current fact nonzero indices (descendantly sorted according to absolute values in L)
				list<pair<int,int>> fact_nz_inds;
				std::function<int(const pair<int, int>& a, const pair<int, int>& b, Faust::MatDense<FPP4,DEVICE>& L_low)> fact_nz_inds_sort_func;

			public:
				void max_L();
			protected:

		/**
		 * After selecting a pivot in iteration previous steps, it's necessary to remove them.
		 * That's what this function is responsible for.
		 * The pivots are nonzero elements of L and their indices are sorted in fact_nz_inds (ascendently according to the pivot values).
		 */
		void update_fact_nz_inds(int p, int q);

		/**
		 * This function is responsible to compute all the coefficients (excepting the identity part)
		 * of the current iteration factor (the coefficients are those from the rotation matrices).
		 */
		void loop_update_fact();

		/**
		 * Chooses the next pivot to zero in diagonalizing process.
		 */
		void choose_pivot();
		/**
		 * Computes the coefficients of the last selected rotation matrix to be put later in current iteration factor.
		 */
		void update_L(Faust::MatSparse<FPP,Cpu> &L);

		/**
		 * Constructs the current factor after computing all the coefficients (of rotation matrices) in temporary buffers (see update_fact()).
		 */
		void finish_fact();

		/**
		 * Constructor.
		 *
		 *
		 * \param t the maximum number of 2D rotation matrices (Givens matrices) to insert in each factor. The effective number can be less than t if there is not enough pivot candidates left in the matrix L of the current iteration.
		 * \param algo instance of GivensFGFT or GivensFGFTComplex (subclasses of GivensFGFTGen).
		 */
		GivensFGFTParallelGen(int t, Faust::GivensFGFTGen<FPP, DEVICE, FPP2, FPP4> & alg);
	};

#include "faust_GivensFGFTParallelGen.hpp"


}
#endif
