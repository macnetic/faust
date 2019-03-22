#include "faust_ParamsPalmFGFT.h"
#include "faust_Palm4MSA.h"
#include "faust_ParamsFGFT.h"

#ifndef __FAUST_PALM4MSA_FGFT_H__
#define __FAUST_PALM4MSA_FGFT_H__

using namespace Faust;

namespace Faust {

	template<typename FPP, Device DEVICE, typename FPP2 = double>
		class Palm4MSAFGFT : public Palm4MSA<FPP, DEVICE, FPP2>
	{
		MatSparse<FPP, DEVICE> D;
		MatDense<FPP,DEVICE> D_grad_over_c;
		public:
			Palm4MSAFGFT(const ParamsPalmFGFT<FPP, DEVICE, FPP2>& params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal=false);
			Palm4MSAFGFT(const MatDense<FPP,DEVICE>& Lap, const ParamsFGFT<FPP,DEVICE,FPP2> & params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal);
			void next_step();
			const MatSparse<FPP, DEVICE>& get_D();
			void get_D(FPP* diag_data);
		private:
			void compute_grad_over_c();
			void compute_lambda();
			void compute_D();
			void compute_D_grad_over_c();
			void compute_c();
	};

#include "faust_Palm4MSAFGFT.hpp"
}


#endif
