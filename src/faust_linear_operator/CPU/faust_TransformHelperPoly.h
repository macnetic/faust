#ifndef __FAUST_TRANSFORM_HELPER_POLY__
#define __FAUST_TRANSFORM_HELPER_POLY__
#include "faust_TransformHelper.h"
namespace Faust
{
	/**
	 * \brief This class aims to represent a Chebyshev polynomial basis as a "Faust".
	 *
	 * It overrides some operations to optimize the performance by taking into account the factors specific structure.
	 */
	template<typename FPP>
		class TransformHelperPoly : public TransformHelper<FPP,Cpu>
		{

			public:
			MatSparse<FPP, Cpu> L; //TODO: must be private (and basisChebyshev a friend of this class)
			MatSparse<FPP, Cpu> twoL; //TODO: must be private
			TransformHelperPoly(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
					const FPP lambda_, const bool optimizedCopy, const bool cloning_fact,
					const bool internal_call) : TransformHelper<FPP,Cpu>(facts, lambda_, optimizedCopy, cloning_fact, internal_call) {}

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x, const bool transpose=false, const bool conjugate=false);
//						MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> A, const bool transpose=false, const bool conjugate=false);
		};

	template<typename FPP>
		TransformHelper<FPP, Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K);
}
#include "faust_TransformHelperPoly.hpp"
#endif

