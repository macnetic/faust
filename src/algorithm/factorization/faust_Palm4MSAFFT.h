#include "faust_ParamsPalmFFT.h"
#include "faust_Palm4MSA.h"

#ifndef __FAUST_PALM4MSA_FFT_H__
#define __FAUST_PALM4MSA_FFT_H__

using namespace Faust;

namespace Faust {

	template<typename FPP, Device DEVICE, typename FPP2 = double>
		class Palm4MSAFFT : public Palm4MSA<FPP, DEVICE, FPP2>
	{
		MatDense<FPP, DEVICE> D; //TODO: later it will need to be Sparse (which needs to add a prototype overload for multiplication in faust_linear_algebra.h)
		Faust::MatDense<FPP,DEVICE> D_grad_over_c; //TODO: move to sparse mat later
		public:
			//TODO: another ctor (like in Palm4MSA) for hierarchical algo. use
			Palm4MSAFFT(const ParamsPalmFFT<FPP, DEVICE, FPP2>& params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal=false);
			virtual void next_step();
		private:
			virtual void compute_grad_over_c();
			virtual void compute_lambda();
			void compute_D();
			void compute_D_grad_over_c();
	};

#include "faust_Palm4MSAFFT.hpp"
}


#endif
