#ifndef __FAUST_TRANSFORM_HELPER_POLY__
#define __FAUST_TRANSFORM_HELPER_POLY__
#include "faust_TransformHelper.h"
namespace Faust
{
	template<typename FPP>
		TransformHelper<FPP, Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0=nullptr, bool lazy_instantiation=true);

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
			MatSparse<FPP, Cpu> *rR;
			std::vector<bool> is_fact_created;
			public:
			TransformHelperPoly(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
					const FPP lambda_, const bool optimizedCopy, const bool cloning_fact,
					const bool internal_call) : TransformHelper<FPP,Cpu>(facts, lambda_, optimizedCopy, cloning_fact, internal_call) {}

			faust_unsigned_int getNbRow() const;
			faust_unsigned_int getNbCol() const;
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &x, const bool transpose=false, const bool conjugate=false);
			Vect<FPP,Cpu> multiply(const FPP* x, const bool transpose=false, const bool conjugate=false);
			void multiply(const FPP* x, FPP* y, const bool transpose=false, const bool conjugate=false);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &X, const bool transpose=false, const bool conjugate=false);
			void multiply(const FPP* X, int n, FPP* out, const bool transpose=false, const bool conjugate=false);
			TransformHelper<FPP, Cpu>* next(int K);
			TransformHelper<FPP, Cpu>* next();
			Vect<FPP, Cpu> poly(MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			MatDense<FPP, Cpu> poly(int n, MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			void poly(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out);
			TransformHelper<FPP, Cpu>* polyFaust(const FPP* coeffs);
			void basisChebyshevT0(MatSparse<FPP,Cpu>* T0=nullptr);
			void basisChebyshevT1();
			void basisChebyshevT2();
			void basisChebyshevTi(int i);
			void basisChebyshev_all();
			void create_rR(const MatSparse<FPP,Cpu>* L);
			~TransformHelperPoly();
			friend TransformHelper<FPP,Cpu>* basisChebyshev<>(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0, bool lazy_instantiation);
		};


}
#include "faust_TransformHelperPoly.hpp"
#endif

