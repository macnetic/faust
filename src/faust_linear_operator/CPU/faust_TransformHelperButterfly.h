#ifndef __FAUST_TRANSFORM_HELPER_DFT__
#define __FAUST_TRANSFORM_HELPER_DFT__
#include "faust_TransformHelper.h"

namespace Faust
{
	template<typename FPP, FDevice DEV>
		class TransformHelperButterfly;

	template<typename FPP>
	class ButterflyMat;

	template<typename FPP>
		class TransformHelperButterfly<FPP, Cpu> : public TransformHelper<FPP, Cpu>
		{
			Vect<FPP, Cpu> perm_d;
			FPP *perm_d_ptr;
			std::vector<unsigned int> bitrev_perm;
			std::vector<ButterflyMat<FPP>> opt_factors;

			// private ctor
			TransformHelperButterfly<FPP, Cpu>(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);


			public:
			static TransformHelperButterfly<FPP,Cpu>* fourierFaust(unsigned int n, const bool norma=true);
			Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu>& x);
			void multiply(const FPP* x, FPP* y);
			Vect<FPP,Cpu> multiply(const FPP* x);
		};

	template<typename FPP>
	class ButterflyMat
	{

		Vect<FPP, Cpu> d1;
		Vect<FPP, Cpu> d2;
		std::vector<int> subdiag_ids;
		int level;

		// \param level: is a 0-base index.
		public:
		ButterflyMat<FPP>(const MatSparse<FPP, Cpu> &factor, int level);

		Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu>& x) const;
		void multiply(const FPP* x, FPP* y, size_t size) const;
		void Display() const;

	};
}
#include "faust_TransformHelperButterfly.hpp"
#endif
