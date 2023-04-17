
#ifndef __GIVENS_FGFT_COMPLEX__
#define __GIVENS_FGFT_COMPLEX__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include "faust_Transform.h"
#include "faust_EigTJGen.h"
#include <cfloat>
#include <vector>

namespace Faust
{

	template<typename FPP, FDevice DEVICE, typename FPP2 = Real<FPP>>
		class EigTJComplex : public EigTJGen<typename FPP::value_type, DEVICE, FPP2, FPP>{
			/**
			 * \class EigTJComplex
			 *
			 * \brief This class implements the Givens FGFT algorithm (Truncated Jacobi algorithm) for the complex matrix case (ideal case being the Hermitian matrix case).
			 * This algorithm is based on the classical Jacobi eigenvalues algorithm.
			 *
			 * See parent class EigTJGen for documentation about members.
			 *
			 *  References:
			 *
			 *  [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
			 *    graph Fourier transforms via multi-layer sparse approximations",
			 *    submitted to IEEE Transactions on Signal and Information Processing
			 *    over Networks.
			 *    <https://hal.inria.fr/hal-01416110>
			 *
			 *    For the complex version of the algorithm, see generalization here:
			 *    https://en.wikipedia.org/wiki/Jacobi_method_for_complex_Hermitian_matrices
			 *
			 */
			/** \brief Temporary storage matrix for maximization of L. */
			MatDense<FPP,DEVICE> C;
			/** \brief Column vector for the rowwise minimization of C (i.e. maximization of L). */
			Vect<FPP,DEVICE> C_min_row;
			public:
				using DType = typename FPP::value_type;
			protected:
				/** \brief Rotation angle for the current iteration's "Givens" matrix. */
				FPP theta1, theta2;

				/**
				 * Function pointer to any step of the algorithm (internal purpose only).
				 */
				typedef void (EigTJComplex<FPP,DEVICE,FPP2>::*substep_fun)();


			public:

				const static unsigned int ERROR_CALC_PERIOD = 100;
				EigTJComplex(MatSparse<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false, const int err_period=100);
				EigTJComplex(MatDense<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false, const int err_period=100);

				virtual void next_step();

				virtual void update_L();

			protected:
				virtual void choose_pivot();

				virtual void max_L();

				void calc_theta();

				/**
				 * This function verifies the sign of theta2 and recompute the Givens coefficients if needed (in case the pivot image is not null). This function is called from update_fact().
				 */
				void check_pivot_image(FPP& c_pp, FPP& c_pq, FPP& c_qp, FPP& c_qq);

				virtual void update_fact();

				virtual void update_L(MatDense<FPP,Cpu> &);

				virtual void update_L(MatSparse<FPP,Cpu> &);

				/**
				 * Computes the first S'*L (only used by update_L() in optimization enabled code).
				 *
				 */
				void update_L_first(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatDense<FPP,DEVICE> & L);

				/**
				 * Computes L*S (only used by update_L() in optimization enabled code).
				 * Must be called after update_L_first() to finish the update of L = S'*L*S.
				 */
				void update_L_second(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatDense<FPP,DEVICE> & L);

				void update_L_first(Eigen::SparseMatrix<FPP, Eigen::RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatSparse<FPP,DEVICE> & L);
				void update_L_second(Eigen::SparseMatrix<FPP, Eigen::RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, Eigen::RowMajor>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatSparse<FPP,DEVICE> & L);

			public:

				const MatSparse<FPP,DEVICE> get_Dspm(const bool ord=false);

		};

}
#include "faust_EigTJComplex.hpp"
#endif
