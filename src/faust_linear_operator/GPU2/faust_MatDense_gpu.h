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
				MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const FPP* data = nullptr);

//				multiply(const Vect<Cpu,FPP> &vec);
				void multiply(MatDense<FPP, GPU2> &other, const char op_this='N');
				void multiply(MatDense<FPP,Cpu> &other, const char op_this='N');
//				void multiply(MatSparse<FPP, Cpu> &other, MatDense<FPP, GPU2>& output, const char op_this='N');
				void multiply(const MatSparse<FPP, Cpu> &other, MatDense<FPP, Cpu>& output, const char op_this='N');
				void multiply(const MatSparse<FPP, Cpu> &other, MatDense<FPP, GPU2>& output, const char op_this='N');
				MatDense<FPP, Cpu> tocpu();
				~MatDense<FPP, GPU2>();
			private:
				static void* dsm_funcs;
				static void* spm_funcs;
				gm_DenseMat_t gpu_mat;
		};

	template <typename FPP>
		void* Faust::MatDense<FPP,GPU2>::dsm_funcs = nullptr;

	template <typename FPP>
		void* Faust::MatDense<FPP,GPU2>::spm_funcs = nullptr;
};
#include "faust_MatDense_gpu_double.hpp"
#endif
#endif
