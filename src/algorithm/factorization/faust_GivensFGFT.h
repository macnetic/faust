
#ifndef __GIVENS_FGFT__
#define __GIVENS_FGFT__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include "faust_Transform.h"
#include "faust_GivensFGFTGen.h"
#include <cfloat>
#include <vector>

namespace Faust
{


	template<typename FPP, FDevice DEVICE, typename FPP2 = Real<FPP>>
		class GivensFGFT : public GivensFGFTGen<FPP, DEVICE, FPP2> {
			/**
			 * \class GivensFGFT
			 *
			 * \brief This class implements the Givens FGFT algorithm.
			 * This algorithm is based on the classical Jacobi eigenvalues algorithm.
			 *
			 * See parent class GivensFGFTGen for documentation about members.
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
			/** \brief Temporary storage matrix for maximization of L. */
			MatDense<FPP,DEVICE> C;
			/** \brief Column vector for the rowwise minimization of C (i.e. maximization of L). */
			Vect<FPP,DEVICE> C_min_row;
			protected:
				/** \brief Rotation angle theta for the current iteration's Givens matrix. */
				FPP2 theta;

				/** \brief In calc_theta() two values are calculated for theta, this boolean is set to true to always choose theta2 (useful for GivensFGFTParallel). */
				bool always_theta2;

				/**
				 * Function pointer to any step of the algorithm (internal purpose only).
				 */
				typedef void (GivensFGFT<FPP,DEVICE,FPP2>::*substep_fun)();


			public:

				const static unsigned int ERROR_CALC_PERIOD = 100;
				/** Algorithm class constructor.
				 * \param Lap The Laplacian matrix to approximate/diagonalize.
				 * \param J The number of iterations, Givens rotations factors.
				 * */
				GivensFGFT(MatSparse<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false);
				GivensFGFT(MatDense<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false);
				/** Destructor */
				virtual ~GivensFGFT() {/*delete[] q_candidates; delete L;*/};

				/**
				 * \brief Algo. main step.
				 *
				 * This function calls all step functions of the algorithm in the proper order to execute only one iteration.
				 *
				 * External code may call this function instead of compute_facts() in order to debug iterations taken one by one.
				 */
				virtual void next_step();

				/**
				 * \brief Algo. step 2.4
				 *
				 * Updates L after Givens factor update for the next iteration.
				 */
				virtual void update_L();

			protected:
				virtual void choose_pivot();

				virtual void max_L();

				void calc_theta();

				virtual void update_fact();

				virtual void update_L(MatDense<FPP,Cpu> &);

				virtual void update_L(MatSparse<FPP,Cpu> &);

				/**
				 * Computes the first S'*L (only used by update_L() in optimization enabled code).
				 *
				 */
				void update_L_first(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatDense<FPP,DEVICE> & L);

				/**
				 * Computes L*S (only used by update_L() in optimization enabled code).
				 * Must be called after update_L_first() to finish the update of L = S'*L*S.
				 */
				void update_L_second(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatDense<FPP,DEVICE> & L);

				void update_L_first(Eigen::SparseMatrix<FPP, Eigen::RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatSparse<FPP,DEVICE> & L);
				void update_L_second(Eigen::SparseMatrix<FPP, Eigen::RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatSparse<FPP,DEVICE> & L);

				void update_err();

			public:


				const MatSparse<FPP,DEVICE> get_Dspm(const bool ord=false);



		};

}
#include "faust_GivensFGFT.hpp"
#endif
