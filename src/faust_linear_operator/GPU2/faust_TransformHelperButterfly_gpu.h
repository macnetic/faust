#ifndef __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#define __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#include "faust_MatButterfly_gpu.h"

namespace Faust
{
	template<typename FPP, FDevice DEV>
		class TransformHelperButterfly;

	template<typename FPP>
		class TransformHelperButterfly<FPP, GPU2> : public TransformHelper<FPP, GPU2>
		{
			int* perm_ids;
			Vect<FPP, GPU2> d_perm;
			bool has_permutation;
			std::vector<MatButterfly<FPP, GPU2>> opt_factors;


			// private ctor
			TransformHelperButterfly(const std::vector<MatGeneric<FPP, Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);

			public:
			~TransformHelperButterfly() { if(perm_ids != nullptr) delete[] perm_ids;}
			static TransformHelper<FPP,GPU2>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,GPU2>* optFaust(const TransformHelper<FPP, GPU2>* F) { throw std::runtime_error("Not yet implemented on GPU");};
			Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu>& x);
			void multiply(const FPP* x, FPP* y);
			Vect<FPP,Cpu> multiply(const FPP* x);
			void multiply(const FPP* A, int A_ncols, FPP* C);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &A);
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A);

		};

}
#include "faust_TransformHelperButterfly_gpu.hpp" 
#endif
#endif
