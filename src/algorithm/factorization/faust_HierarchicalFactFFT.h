#ifndef __FAUST_HIERARCHICAL_FACT_FFT_H__
#define __FAUST_HIERARCHICAL_FACT_FFT_H__
#include "faust_Palm4MSAFFT.h"
#include "faust_ParamsFFT.h"
#include "faust_HierarchicalFact.h"

using namespace Faust;

namespace Faust
{
    template<typename FPP,Device DEVICE,typename FPP2 = double>
    class HierarchicalFactFFT : public HierarchicalFact<FPP, DEVICE, FPP2>
    {

		Palm4MSAFFT<FPP,DEVICE,FPP2> palm_global;
		static const char * m_className;

		public:
		//TODO: move def. code in .hpp
		HierarchicalFactFFT(const MatDense<FPP,DEVICE>& U, const MatDense<FPP,DEVICE>& Lap,const ParamsFFT<FPP,DEVICE,FPP2>& params, BlasHandle<DEVICE> cublasHandle, SpBlasHandle<DEVICE> cusparseHandle): HierarchicalFact<FPP, DEVICE, FPP2>(U, params, cublasHandle, cusparseHandle), palm_global(Palm4MSAFFT<FPP,DEVICE,FPP2>(Lap, params, cublasHandle, true))
																																			//TODO: verify if palm_global is really initialized after parent ctor call
		{
			if ((U.getNbRow() != params.m_nbRow) |  (U.getNbCol() != params.m_nbCol))
				handleError(m_className,"constructor : params and Fourier matrix U haven't compatible size");
			if((Lap.getNbRow() != params.m_nbRow) |  (Lap.getNbCol() != params.m_nbCol))
				handleError(m_className,"constructor : params and Laplacian matrix Lap haven't compatible size");
		}


	};
}

#include "faust_HierarchicalFactFFT.hpp"

#endif
