#ifndef __FAUST_MATDENSE_GPU2__
#define __FAUST_MATDENSE_GPU2__
#ifdef USE_GPU_MOD
#include "faust_MatDense.h"
#define __GM_LOADER__
#include "gm_interf.h"
namespace Faust
{

	template<typename FPP,FDevice DEVICE>
    class MatDense;

	template<typename FPP>
		class MatDense<FPP, GPU2> : MatDense<FPP, Cpu>
		{
			public:
				MatDense(const FPP* data, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol);

//				multiply(const MatDense<Cpu,FPP> &other);
//				multiply(const MatSparse<Cpu,FPP> &other);
//				multiply(const Vect<Cpu,FPP> &vec);
				void multiply(MatDense<FPP, GPU2> &other, const char op_this='N');
				MatDense<FPP, Cpu> tocpu();
			private:
				static void* dsm_funcs;
				gm_DenseMat_t gpu_mat;
		};

	template <typename FPP>
		void* Faust::MatDense<FPP,GPU2>::dsm_funcs = nullptr;

};
#include "faust_MatDense_gpu_double.hpp"
#endif
#endif
