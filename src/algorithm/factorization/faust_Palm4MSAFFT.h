#include "faust_ParamsPalmFFT.h"
#include "faust_Palm4MSA.h"

#ifndef __FAUST_PALM4MSA_FFT_H__
#define __FAUST_PALM4MSA_FFT_H__

using namespace Faust;

namespace Faust {

	template<typename FPP, Device DEVICE, typename FPP2 = double>
		class Palm4MSAFFT : public Palm4MSA<FPP, DEVICE, FPP2>
	{
		public:
			//TODO: another ctor (like in Palm4MSA) for hierarchical algo. use
			Palm4MSAFFT(const ParamsPalmFFT<FPP, DEVICE, FPP2>& params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal=false);
		private:
			virtual void compute_grad_over_c();
			virtual void compute_lambda();
	};

#include "faust_Palm4MSAFFT.hpp"
}


#endif
