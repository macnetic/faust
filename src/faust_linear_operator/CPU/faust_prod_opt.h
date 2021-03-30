#ifndef __FAUST_PROD_OPT__
#define __FAUST_PROD_OPT__
#include <vector>
#include "faust_MatDense.h"
#ifdef __APPLE__
#ifdef _MUL_OMP_
#include "omp.h"
#endif
#endif
namespace Faust
{
	enum FaustMulMode
	{
		DEFAULT,
		ORDER_ALL_ENDS,
		ORDER_1ST_BEST,
		ORDER_ALL_BEST_CONVDENSE,
		ORDER_ALL_BEST_MIXED,
		CPP_PROD_PAR_REDUC,
		OMP_PROD_PAR_REDUC,
		TORCH_CPU,
		TORCH_CPU_BEST_ORDER,
		TORCH_CPU_DENSE_ROW_TO_TORCH,
	};

	/**
	 * Computes facts product by varying the chosen end side to multiply two matrices at each stage.
	 * The choice is made for efficiency by choosing the less costly side.
	 *
	 * For example: if facts = {A B C D}, the first step is to choose to compute AB or CD, if the number of scalar product is smaller for AB then AB is computed, otherwise this is CD. For next iteration the same question is considered whatever is the previous choice.
	 * If it was [AB] C D, then the costs of [AB]C and CD are compared. Likewise, if the first calculation was CD, the cost of AB and B[CD] are compared in order to proceed the calculation in the less costly manner.
	 *
	 * The function computes out = alpha*op_0(facts[0])*...*op_n(facts[n])+beta*out.
	 * op_i(facts[i]) is defined by the transconj_flags[i] if i < transconj_flags.size(): 'N' means op_i is nop (non-operation or identity), 'T' op_i is the transpose operation, and 'H' for the adjoint op.
	 * if i >= transconj_flags.size() and i > 1 then op_i is defined by the transconj_flags[j] with j == transconj_flags.size()-1.
	 * It's handy to define the same op for all facts[i], as it's the case with the default argument transconj_flags (all factors are multiplied after a non-op.).
	 */
	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_ends(std::vector<MatDense<FPP,DEVICE>*>& facts, MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	/**
	 *
	 * This function does the same as multiply_order_opt_all_ends but is capable to multiply not only on the ends of the matrix chain but also in the middle if a better complexity guides to this.
	 *
	 */
	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_best(std::vector<MatDense<FPP,DEVICE>*>& facts, MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	/**
	 * This function does the same as multiply_order_opt_all_best except that here facts matrices can be MatDense or MatSparse too. In the latter instead of dimensions this is the nnz that is taken into account to minimize the cost.
	 */
	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_best(std::vector<MatGeneric<FPP,DEVICE>*>& facts, MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	/**
	 * This functions does the same as multiply_order_opt_all_best but contrary to the latter it optimizes the order only for the first multiplication (wherever it is located in the chain matrix). Then it multiplies all matrices on the left and on the right of the chosen "best" position.
	 *
	 * For example: to multiply A B C D E, if BC is the less costly choice to multiply.
	 * It multiplies CD, then A(B(CD)), then (A(B(CD)))E.
	 */
	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_first_best(std::vector<MatDense<FPP,DEVICE>*>& facts, MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	/**
	 * This function is a wrapper to other order based optimization methods/functions (declared above).
	 *
	 */
	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt(const int mode, std::vector<MatGeneric<FPP,DEVICE>*>& facts, MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));


	/**
	 * Computes the data matrices product by proceeding to a parallel reduction.
	 *
	 * Example of the principle: if data.size() == 4 and two threads are available, in the first stage the thread #1 multiplies the matrices 0 and 1, while the thread #2 multiplies the matrices 2 and 3. It results into two new matrices 0_1 and 1_2 that it remains to multiply to obtain the whole product. The thread #1 handles it and that's it.
	 *
	 * NOTE: this method was for test purpose and it turns out that when Eigen multithread support is enabled there is not point to use this parallel reduction from a performance point of view (the parallelization based on vector dot product when multiplying two matrices is quicker than this parallel reduction).
	 */
	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> multiply_omp(const std::vector<MatGeneric<FPP,DEVICE>*> data, const MatDense<FPP,DEVICE> A, const char opThis);

	/**
	 * This function is the run method of thread created by multiply_par declared below.
	 */
	template<typename FPP, FDevice DEVICE>
	void multiply_par_run(int nth, int thid, int num_per_th, int data_size, char opThis, std::vector<MatGeneric<FPP, DEVICE>*>& data, std::vector<MatDense<FPP,DEVICE>*>& mats, std::vector<std::thread*>& threads, std::mutex &, std::condition_variable &, int &);

	/**
	 * This function is the same as multiply_omp except that it is implemented using C++ threads instead of OpenMP.
	 * For the same reason as for multiply_omp this method isn't useful (see NOTE in multiply_omp code doc.).
	 */
	template<typename FPP, FDevice DEVICE>
	MatDense<FPP,DEVICE> multiply_par(const std::vector<MatGeneric<FPP,DEVICE>*>& data, const MatDense<FPP,DEVICE> A, const char opThis);

}
#include "faust_prod_opt.hpp"
#endif
