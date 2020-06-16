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
		PAR_SUBPRODS_CPP,
		PAR_SUBPRODS_OMP,
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

	template<typename FPP, FDevice DEVICE>
		Faust::MatDense<FPP,DEVICE> multiply_omp(const std::vector<Faust::MatGeneric<FPP,DEVICE>*> data, const Faust::MatDense<FPP,DEVICE> A, const char opThis);

	template<typename FPP, FDevice DEVICE>
	void multiply_par_run(int nth, int thid, int num_per_th, int data_size, char opThis, std::vector<Faust::MatGeneric<FPP, DEVICE>*>& data, std::vector<Faust::MatDense<FPP,DEVICE>*>& mats, std::vector<std::thread*>& threads, std::mutex &, std::condition_variable &, int &);

	template<typename FPP, FDevice DEVICE>
	MatDense<FPP,DEVICE> multiply_par(const std::vector<Faust::MatGeneric<FPP,DEVICE>*>& data, const Faust::MatDense<FPP,DEVICE> A, const char opThis);

}
#include "faust_prod_opt.hpp"
#endif
