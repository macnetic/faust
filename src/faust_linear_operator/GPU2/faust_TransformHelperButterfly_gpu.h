#ifndef __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#define __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"

namespace Faust
{
	template<typename FPP, FDevice DEV>
		class TransformHelperButterfly;

	template<typename FPP, FDevice DEV>
	class ButterflyMat;

	template<typename FPP>
		class TransformHelperButterfly<FPP, GPU2> : public TransformHelper<FPP, GPU2>
		{
			int* perm_ids;
			Vect<FPP, GPU2> d_perm;
			bool has_permutation;
			std::vector<ButterflyMat<FPP, GPU2>> opt_factors;


			// private ctor
			TransformHelperButterfly(const std::vector<MatGeneric<FPP, Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);

			~TransformHelperButterfly() { delete[] perm_ids;}
			public:
			static TransformHelper<FPP,GPU2>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,GPU2>* optFaust(const TransformHelper<FPP, GPU2>* F) { throw std::runtime_error("Not yet implemented on GPU");};
			Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu>& x);
			void multiply(const FPP* x, FPP* y);
			Vect<FPP,Cpu> multiply(const FPP* x);
			void multiply(const FPP* A, int A_ncols, FPP* C);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &A);
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A);

		};

	template<typename FPP>
	class ButterflyMat<FPP, GPU2>
	{

		Vect<FPP, GPU2> d1;
		Vect<FPP, GPU2> d2;
		int* subdiag_ids;
#ifdef USE_PYTHONIC
		long *subdiag_ids_ptr;
#endif
		int level;

		// \param level: is a 0-base index.
		public:
		ButterflyMat(const MatSparse<FPP, Cpu> &factor, int level);

		//TODO: constness of multiply member functions
		MatDense<FPP, GPU2> multiply(const FPP* x);
		void Display() const;

		MatDense<FPP, GPU2> multiply(const FPP* A, int A_ncols);
		MatDense<FPP, GPU2> multiply(MatDense<FPP,GPU2> &A);
		void multiply(MatDense<FPP,GPU2> &A, MatDense<FPP, Cpu> & out);
//		MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A);
		const Vect<FPP, GPU2>& getD1() {return d1;};
		const Vect<FPP, GPU2>& getD2() {return d2;};

		~ButterflyMat() { delete[] subdiag_ids;}
	};
}
#include "faust_TransformHelperButterfly_gpu.hpp" //TODO
#endif
#endif
