#ifndef __FAUST_TRANSFORM_HELPER_POLY__
#define __FAUST_TRANSFORM_HELPER_POLY__
#include "faust_TransformHelper.h"
namespace Faust
{
	template<typename FPP>
		TransformHelper<FPP, Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K);

	template<typename FPP>
		void poly(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out);

	/**
	 * \brief This class aims to represent a Chebyshev polynomial basis as a "Faust".
	 *
	 * It overrides some operations to optimize the performance by taking into account the factors specific structure.
	 */
	template<typename FPP>
		class TransformHelperPoly : public TransformHelper<FPP,Cpu>
		{

			static RefManager ref_man;
			MatSparse<FPP, Cpu> *L;
			MatSparse<FPP, Cpu> *twoL;
			MatSparse<FPP, Cpu> *rR;
			public:
			TransformHelperPoly(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
					const FPP lambda_, const bool optimizedCopy, const bool cloning_fact,
					const bool internal_call) : TransformHelper<FPP,Cpu>(facts, lambda_, optimizedCopy, cloning_fact, internal_call) {}

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &x, const bool transpose=false, const bool conjugate=false);
			Vect<FPP,Cpu> multiply(const FPP* x, const bool transpose=false, const bool conjugate=false);
			void multiply(const FPP* x, FPP* y, const bool transpose=false, const bool conjugate=false);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &X, const bool transpose=false, const bool conjugate=false);
			TransformHelper<FPP, Cpu>* next(int K);
			TransformHelper<FPP, Cpu>* next();
			Vect<FPP, Cpu> poly(MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			MatDense<FPP, Cpu> poly(int n, MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			void poly(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out);
			TransformHelper<FPP, Cpu>* polyFaust(const FPP* coeffs);
			~TransformHelperPoly();
			friend TransformHelper<FPP,Cpu>* basisChebyshev<>(MatSparse<FPP,Cpu>* L, int32_t K);
		};


}
#include "faust_TransformHelperPoly.hpp"
#endif

