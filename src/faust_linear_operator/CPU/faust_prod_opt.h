#ifndef __FAUST_PROD_OPT__
#define __FAUST_PROD_OPT__
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
		GPU_MOD
	};

	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_ends(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_all_best(std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt_first_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

	template<typename FPP, FDevice DEVICE>
		void multiply_order_opt(const int mode, std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));


	/**
	 * Computes the data matrices product by proceeding to a parallel reduction.
	 *
	 * Example of the principle: if data.size() == 4 and two threads are available, in the first stage the thread #1 multiplies the matrices 0 and 1, while the thread #2 multiplies the matrices 2 and 3. It results into two new matrices 0_1 and 1_2 that it remains to multiply to obtain the whole product. The thread #1 handles it and that's it. 
	 *
	 * NOTE: this method was for test purpose and it turns out that when Eigen multithread support is enabled there is not point to use this parallel reduction from a performance point of view (the parallelization based on vector dot product when multiplying two matrices is quicker than this parallel reduction).
	 */
	template<typename FPP, FDevice DEVICE>
		Faust::MatDense<FPP,DEVICE> multiply_omp(const std::vector<Faust::MatGeneric<FPP,DEVICE>*> data, const Faust::MatDense<FPP,DEVICE> A, const char opThis);

	/**
	 * This function is the run method of thread created by multiply_par declared below.
	 */
	template<typename FPP, FDevice DEVICE>
	void multiply_par_run(int nth, int thid, int num_per_th, int data_size, char opThis, std::vector<Faust::MatGeneric<FPP, DEVICE>*>& data, std::vector<Faust::MatDense<FPP,DEVICE>*>& mats, std::vector<std::thread*>& threads, std::mutex &, std::condition_variable &, int &);

	/**
	 * This function is the same as multiply_omp except that it is implemented using C++ threads instead of OpenMP.
	 * For the same reason as for multiply_omp this method isn't useful (see NOTE in multiply_omp code doc.).
	 */
	template<typename FPP, FDevice DEVICE>
	MatDense<FPP,DEVICE> multiply_par(const std::vector<Faust::MatGeneric<FPP,DEVICE>*>& data, const Faust::MatDense<FPP,DEVICE> A, const char opThis);

}
#include "faust_prod_opt.hpp"
#endif
