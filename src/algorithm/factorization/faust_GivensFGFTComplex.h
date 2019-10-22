
#ifndef __GIVENS_FGFT__
#define __GIVENS_FGFT__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include "faust_Transform.h"
#include <vector>

namespace Faust {


	template<typename FPP, Device DEVICE, typename FPP2 = float>
		class GivensFGFTComplex {
			/**
			 * \class Faust::GivensFGFTComplex
			 *
			 * \brief This class implements the Givens FGFT algorithm.
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
			/** \brief Temporary storage matrix for maximization of L. */
			Faust::MatDense<FPP,DEVICE> C;
			/** \brief Column vector for the rowwise minimization of C (i.e. maximization of L). */
			Faust::Vect<FPP,DEVICE> C_min_row;
			/** \brief Pivot candidates q coordinates. */
			int* q_candidates;  /* default IndexType for underlying eigen matrix is int. */
			protected:
				const static unsigned int ERROR_CALC_PERIOD = 100;
				/** \brief Fourier matrix factorization matrices (Givens matrix). */
				vector<Faust::MatSparse<FPP,DEVICE>> facts;
				/** \brief Diagonalization approximate of Laplacian. */
				Faust::Vect<FPP,DEVICE> D;
				/** \brief Queue of errors (cf. calc_err()). */
				vector<FPP2> errs;
				/** \brief Pivot choices (p, q) coordinates. */
				vector<pair<int,int>> coord_choices;
				/** \brief Graph Laplacian to diagonalize/approximate. */
				Faust::MatGeneric<FPP, DEVICE>& Lap;
				/** \brief L iteration factor:  L_i = S^T L_{i-1} S, initialized from Lap (with S being facts[i]). */
				Faust::MatGeneric<FPP, DEVICE> *L;
				/** \brief Rotation angle theta for the current iteration's Givens matrix. */
				FPP theta1, theta2;

				/* Precomputed model identity matrix to init. facts[ite] before update.
				 * Identity matrix is completed later with cos/sin coefficients (in update_fact()).
				 */

				/** \brief Defines the rows of facts[ite]. */
				vector<int> fact_mod_row_ids;
				/** \brief Defines the columns of facts[ite]. */
				vector<int> fact_mod_col_ids;
				/** \brief Defines the coefficients of facts[ite]. */
				vector<FPP> fact_mod_values;


				/** \brief Ordered indices of D to get increasing eigenvalues along the diagonal. */
				vector<int> ord_indices;
				/** \brief Cache for the ordered D. */
				Faust::Vect<FPP,DEVICE> ordered_D;
				/** \brief true if D has already been ordered (order_D() was called). */
				bool is_D_ordered;
				int D_order_dir;

				/** \brief The level of verbosity (0 for nothing, 1 for iteration numbers,...) */
				unsigned int verbosity;

				/**
				 * \brief Row index for the selected pivot in L.
				 */
				int p,
					/**
					 * \brief Column index for the selected pivot in L.
					 */q;

					/**
					 * \brief Current iteration number (and index to the current factor to update).
					 */
					unsigned int ite;

				/** \brief In calc_theta() two values are calculated for theta, this boolean is set to true to always choose theta2 (useful for GivensFGFTComplexParallel). */
				bool always_theta2;

				bool stoppingCritIsError;
				double stoppingError;
				bool errIsRel;

				/**
				 * Function pointer to any step of the algorithm (internal purpose only).
				 */
				typedef void (GivensFGFTComplex<FPP,DEVICE,FPP2>::*substep_fun)();


			public:

				/** Algorithm class constructor.
				 * \param Lap The Laplacian matrix to approximate/diagonalize.
				 * \param J The number of iterations, Givens rotations factors.
				 * */
				GivensFGFTComplex(Faust::MatSparse<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true);
				GivensFGFTComplex(Faust::MatDense<FPP,DEVICE>& Lap, int J, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true);
				/** Destructor */
				virtual ~GivensFGFTComplex() {delete[] q_candidates; delete L;};

				/**
				 * \brief Algo. main step.
				 *
				 * This function calls all step functions of the algorithm in the proper order to execute only one iteration.
				 *
				 * External code may call this function instead of compute_facts() in order to debug iterations taken one by one.
				 */
				virtual void next_step();

				/**
				 * \brief Algo. main loop (facts.size() iterations).
				 *
				 * It's the entry point to user code.
				 * It calls next_step() in order to execute one iteration.
				 *
				 */
				void compute_facts();

				/**
				 * \brief Algo. step 2.4
				 *
				 * Updates L after Givens factor update for the next iteration.
				 */
				virtual void update_L();

			protected:
				/** \brief Algo. step 2.1.
				*/
				virtual void choose_pivot();

				/** \brief Algo. step 2.1.1
				 *
				 *	Computes the max of L or sorts it.
				 */
				virtual void max_L();

				/** \brief Algo. step 2.2.
				 *
				 * Computes theta angle for the current Givens factor.
				 *
				 */
				void calc_theta();

				/**
				 * \brief Algo. step 2.3.
				 *
				 * Updates the current Givens factor according to the pivot chosen for the iteration.
				 *
				 */
				virtual void update_fact();

				virtual void update_L(MatDense<FPP,Cpu> &);

				virtual void update_L(MatSparse<FPP,Cpu> &);

				/**
				 * Computes the first S'*L (only used by update_L() in optimization enabled code).
				 *
				 */
				void update_L_first(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatDense<FPP,DEVICE> & L);

				/**
				 * Computes L*S (only used by update_L() in optimization enabled code).
				 * Must be called after update_L_first() to finish the update of L = S'*L*S.
				 */
				void update_L_second(Faust::Vect<FPP,DEVICE>& L_vec_p, Faust::Vect<FPP,DEVICE>& L_vec_q, const FPP& c_pp, const FPP& c_pq, const FPP& c_qp, const FPP& c_qq, int p, int q, MatDense<FPP,DEVICE> & L);

				void update_L_first(Eigen::SparseMatrix<FPP, RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP& c, const FPP& s, int p, int q, Faust::MatSparse<FPP,DEVICE> & L);
				void update_L_second(Eigen::SparseMatrix<FPP, RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP& c, const FPP& s, int p, int q, Faust::MatSparse<FPP,DEVICE> & L);
				/**
				 * \brief Algo. step 2.5.
				 *
				 * Updates the diagonal approximate D from L (the diagonal part is taken).
				 *
				 */
				void update_D();

				/**
				 * \brief Algo. step 2.6.
				 *
				 * Computes the error of approximation for the current iteration.
				 *
				 */
				void update_err();

				/**
				 * Order (sort) D ascendantly into ordered_D and keeps ordered indices in ord_indices.
				 */
				void order_D();

				/**
				 * Order D into ascendantly or descendantly into order_D and keeps ordered indices in ord_indices.
				 *
				 * \param order -1 to order descendantly, 1 to order ascendantly.
				 */
				void order_D(int order);


			public:

				/**
				 * Returns the ordered indices of D to get increasing eigenvalues along the diagonal.
				 *
				 * @see order_D()
				 */
				const vector<int>& get_ord_indices();


				/**
				 * Returns the vector of errors computed by calc_err() during the algorithm iterations.
				 */
				const vector<FPP2>& get_errs() const;

				/**
				 * Returns the specific j-th iteration's error (computed in calc_err()).
				 */
				FPP2 get_err(int j) const;

				/**
				 * Returns the diag. vector D in its current status (which is updated at each iteration).
				 */
				const Faust::Vect<FPP,DEVICE>& get_D(const bool ord=false);

				/**
				 * \param ord 0 (no sort), 1 (ascendant sort), -1 (descendant sort).
				 */
				const Faust::Vect<FPP,DEVICE>& get_D(const int ord=0);


				/** \brief Returns the diagonal vector as a sparse matrix.
				 *
				 * \note This function makes copies, is not intented for repeated use (use get_D() for optimized calls).
				 *
				 **/
				const Faust::MatSparse<FPP,DEVICE> get_Dspm(const bool ord=false);

				/**
				 * Returns the diagonal by copying it in the buffer diag_data (should be allocated by the callee).
				 *
				 * \param ord true to get the data ordered by ascendant eigenvalues false otherwise.
				 */
				void get_D(FPP* diag_data, const bool ord=false);

				/**
				 *
				 * \param ord: 0 for no sort, -1 for descending order, 1 for descending order.
				 */
				void get_D(FPP* diag_data, const int ord=0);


				/**
				 * Computes and returns the Fourier matrix approximate from the Givens factors computed up to this time.
				 *
				 */
				const MatDense<FPP,DEVICE> compute_fourier(const bool ord=false);

				/**
				 * Returns the matrix L.
				 *
				 * @note this matrix is not the Laplacian matrix from which the algorithm has started from. This is a matrix initialized to the Laplacian before the first iteration and then is computed such that L = S^t L. S being the current Givens factor matrix updated (i.e. facts[ite]).
				 *
				 * @see get_Lap()
				 *
				 */
				const MatGeneric<FPP,DEVICE>& get_L() const ;

				/**
				 * Returns the vector of all pivot choices (p, q) coordinates made by the algorithm until the last iteration.
				 */
				const vector<pair<int,int>>& get_coord_choices() const;

				/**
				 * Returns the j-th iteration's pivot choice (p, q).
				 */
				void get_coord_choice(int j, int& p, int& q) const;

				/**
				 * Returns the Laplacian matrix unmodified (as it was when passed to the constructor).
				 */
				const MatDense<FPP,DEVICE>& get_Lap() const;

				/**
				 * Returns the vector of Givens matrices at this stage of algorithm execution (terminated or not).
				 */
				const vector<MatSparse<FPP,DEVICE>>& get_facts() const;


				/**
				 * Returns a Faust::Transform object with copy of facts into it.
				 *
				 * \param ord true to get the Transform's facts ordering the last one columns according to ascendant eigenvalues, false to let facts as they went out from the algorithm (without reordering).
				 */
				Faust::Transform<FPP,DEVICE> get_transform(bool ord);
				/**
				 * \param ord -1 for descending order, 1 for ascending order.
				 */
				Faust::Transform<FPP,DEVICE> get_transform(int ord);


		};

}
#include "faust_GivensFGFTComplex.hpp"
#endif
