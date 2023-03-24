
#ifndef __GIVENS_FGFT_GEN__
#define __GIVENS_FGFT_GEN__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include "faust_Transform.h"
#include <cfloat>
#include <vector>

namespace Faust
{

	template<typename FPP, FDevice DEVICE, typename FPP2, typename FPP4>
		class GivensFGFTParallelGen;

	template<typename FPP, FDevice DEVICE, typename FPP2 = Real<FPP>, typename FPP4 = FPP>
		class GivensFGFTGen {
			/**
			 * \class GivensFGFTGen
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
			friend class GivensFGFTParallelGen<FPP, DEVICE, FPP2, FPP4>;
			/** \brief Temporary storage matrix for maximization of L. */
			//			MatDense<FPP4,DEVICE> C;
			/** \brief Column vector for the rowwise minimization of C (i.e. maximization of L). */
			//			Vect<FPP4,DEVICE> C_min_row;
			protected:
			/** \brief Fourier matrix/eigenvectors factorization matrices (Givens matrices). */
			std::vector<MatSparse<FPP4,DEVICE>> facts;

			/** \brief L iteration factor:  L_i = S^T L_{i-1} S, initialized from Lap (with S being facts[i]). */
			MatGeneric<FPP4, DEVICE>* L;
			/** \brief Pivot candidates q coordinates. */
			int* q_candidates;  /* default IndexType for underlying eigen matrix is int. */

			/** \brief The number of targeted transform factors (when only one Givens rotation is stored per factor).*/
			int J;
			/** \brief approximate eigenvalues vector. */
			Vect<FPP,DEVICE> D;
			/** \brief Queue of errors (cf. update_err()). */
			vector<FPP2> errs;
			/** \brief Pivot choices (p, q) coordinates. */
			vector<pair<int,int>> coord_choices;
			/** \brief Graph Laplacian to diagonalize/approximate. */
			MatGeneric<FPP4, DEVICE>& Lap;
			/** \brief Laplacian dimension */
			unsigned int dim_size;
			/** \brief Laplacian Frobenius norm */
			FPP2 Lap_squared_fro_norm;
			/** \brief Rotation angle theta for the current iteration's Givens matrix. */
			//				FPP2 theta;

			/* Precomputed model identity matrix to init. facts[ite] before update.
			 * Identity matrix is completed later with cos/sin coefficients (in update_fact()).
			 */

			/** \brief Defines the rows of facts[ite]. */
			vector<int> fact_mod_row_ids;
			/** \brief Defines the columns of facts[ite]. */
			vector<int> fact_mod_col_ids;
			/** \brief Defines the coefficients of facts[ite]. */
			vector<FPP4> fact_mod_values;


			/** \brief Ordered indices of D to get increasing eigenvalues along the diagonal. */
			vector<int> ord_indices;

			/** \brief inverse permutation of ord_indices (needed to retrieve start undefined order). */
			vector<int> inv_ord_indices;
			/** \brief Cache for the ordered D. */
			Vect<FPP,DEVICE> ordered_D;
			/** \brief true if D has already been ordered (order_D() was called). */
			bool is_D_ordered;
			/** \brief 1 if eigenvalues has been ordered in ascending order -1 ortherwise (other values is for undefined order). */
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

			/** \brief true if the stopping criterion is error (otherwise it's the number of iterations/number of Givens matrices) */
			bool stoppingCritIsError;
			/** \brief error value according to algorithm stops if stoppingCritIsError == true */
			double stoppingError;
			/** \brief true if the stopping error is taken as relative error (absolute otherwise). */
			bool errIsRel;
			/** \brief (false to default) true to force the computation of the transform even if it doesn't worth it in term of complexity */
			bool enable_large_Faust;
			/** \brief Period (in number of iterations) according to the error is calculated (it applies if stoppingCritIsError == true or verbosity > 2)
			*/
			int err_period;

			/** \brief An arbitrary tag to identify more particularly this instance of the algorithm in output */
			std::string tag;

			public:

				/** Algorithm class constructor.
				 * \param Lap The Laplacian matrix to approximate/diagonalize.
				 * \param J The number of iterations, Givens rotations factors.
				 * TODO: complete argument list
				 * */
				// the MatSparse/MatDense constructor rely on this one
				GivensFGFTGen(MatGeneric<FPP4,DEVICE>* Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust = false, const int err_period=100);

				GivensFGFTGen(MatSparse<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust = false, const int err_period=100);
				GivensFGFTGen(MatDense<FPP4, DEVICE> & Lap, int J, unsigned int verbosity /* deft val == 0 */, const double stoppingError, const bool errIsRel, const bool enable_large_Faust = false, const int err_period=100);

				/** Destructor */
				virtual ~GivensFGFTGen() {delete[] q_candidates; delete L;};

				/**
				 * \brief Algo. main step.
				 *
				 * This function calls all step functions of the algorithm in the proper order to execute only one iteration.
				 *
				 * External code may call this function instead of compute_facts() in order to debug iterations taken one by one.
				 */
				virtual void next_step()=0;

				/**
				 * \brief Algo. main loop (facts.size() iterations).
				 *
				 * It's the entry point to user code.
				 * It calls next_step() in order to execute one iteration.
				 *
				 */
				virtual void compute_facts();

				/**
				 * \brief Algo. step 2.4
				 *
				 * Updates L after Givens factor update for the next iteration.
				 */
				virtual void update_L()=0;

			protected:
				/** \brief Algo. step 2.1.
				*/
				virtual void choose_pivot()=0;

				/** \brief Algo. step 2.1.1
				 *
				 *	Computes the max of L or sorts it.
				 */
				virtual void max_L()=0;

				/** \brief Algo. step 2.2.
				 *
				 * Computes theta angle for the current Givens factor.
				 *
				 */
				virtual void calc_theta()=0;

				/**
				 * \brief Algo. step 2.3.
				 *
				 * Updates the current Givens factor according to the pivot chosen for the iteration.
				 *
				 */
				virtual void update_fact()=0;

				//				virtual void update_L(MatDense<FPP,Cpu> &);
				//
				//				virtual void update_L(MatSparse<FPP,Cpu> &);

				/**
				 * Computes the first S'*L (only used by update_L() in optimization enabled code).
				 *
				 */
				//				void update_L_first(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatDense<FPP,DEVICE> & L);

				/**
				 * Computes L*S (only used by update_L() in optimization enabled code).
				 * Must be called after update_L_first() to finish the update of L = S'*L*S.
				 */
				//				void update_L_second(Vect<FPP,DEVICE>& L_vec_p, Vect<FPP,DEVICE>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatDense<FPP,DEVICE> & L);

				//				void update_L_first(Eigen::SparseMatrix<FPP, RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatSparse<FPP,DEVICE> & L);
				//				void update_L_second(Eigen::SparseMatrix<FPP, RowMajor> & L_vec_p, Eigen::SparseMatrix<FPP, RowMajor>& L_vec_q, const FPP2& c, const FPP2& s, int p, int q, MatSparse<FPP,DEVICE> & L);
				/**
				 * \brief Algo. step 2.5.
				 *
				 * Updates the diagonal approximate D from L (the diagonal part is taken).
				 *
				 */
				virtual void update_D();

				/**
				 * \brief Algo. step 2.6.
				 *
				 * Computes the error of approximation for the current iteration.
				 *
				 */
				virtual void update_err();

				/**
				 * Sort D in descending order into ordered_D and keeps ordered indices in ord_indices.
				 */
				void order_D();

				/**
				 * Sort D into order_D and keeps ordered indices in ord_indices.
				 *
				 * \param order -1 for descending order, 1 to ascending order.
				 */
				void order_D(int order);


			public:

				/**
				 *  Computes the error norm of D compared to L.
				 *  It can be relative or absolute depending on this->errIsRel.
				 */
				FPP2 calc_err();

				/**
				 * Returns the ordered indices of D to get increasing eigenvalues along the diagonal.
				 *
				 * @see order_D()
				 */
				const vector<int>& get_ord_indices(const int order=1);


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
				const Vect<FPP,DEVICE>& get_D(const bool ord=false);

				/**
				 * \param ord 0 (no sort), 1 (ascendant sort), -1 (descendant sort).
				 */
				const Vect<FPP,DEVICE>& get_D(const int ord=0);

				/** \brief Returns the diagonal vector as a sparse matrix.
				 *
				 * \note This function makes copies, is not intented for repeated use (use get_D() for optimized calls).
				 *
				 **/
				template<typename FPP3>
					void get_Dspm(MatSparse<FPP3,DEVICE>&, const bool ord=false);

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
				const MatDense<FPP4,DEVICE> compute_fourier(const bool ord=false);

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
				const MatDense<FPP4,DEVICE>& get_Lap() const;

				/**
				 * Returns the vector of Givens matrices at this stage of algorithm execution (terminated or not).
				 */
				const vector<MatSparse<FPP4,DEVICE>>& get_facts() const;

				/**
				 * Number of Givens factors used for the eigenvector approximate so far.
				 */
				const size_t nfacts() const
				{
					return this->get_facts().size();
				}

				void set_tag(const std::string &tag)
				{
					this->tag = tag;
				}

				std::string get_tag() const
				{
					return tag;
				}


				/**
				 * Returns a Transform object with copy of facts into it.
				 *
				 * \param ord true to get the Transform's facts ordering the last one columns according to ascending eigenvalues, false to let facts as they went out from the algorithm (without reordering).
				 */
				Transform<FPP4,DEVICE> get_transform(const bool ord);
				/**
				 * Returns the eigen vectors as a Transform t of Givens matrices.
				 *
				 * \param ord -1 for descending order, 1 for ascending order.
				 *
				 * Warning: if copy == false and ord == true the caller is responsible to free the last permuted factor of t.
				 */
				Transform<FPP4,DEVICE> get_transform(const int ord, const bool copy=true, const int first_nfacts=-1);


		};

}
#include "faust_GivensFGFTGen.hpp"
#endif
