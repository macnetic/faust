#ifndef __FAUST_HIERARCHICAL_FACT_FGFT_H__
#define __FAUST_HIERARCHICAL_FACT_FGFT_H__
#include "faust_Palm4MSAFGFT.h"
#include "faust_ParamsFGFT.h"
#include "faust_HierarchicalFact.h"

using namespace Faust;

namespace Faust
{
    template<typename FPP,FDevice DEVICE,typename FPP2 = double>
    class HierarchicalFactFGFT : public HierarchicalFact<FPP, DEVICE, FPP2>
    {

		static const char * m_className;

		public:
		//TODO: move def. code in .hpp
		HierarchicalFactFGFT(const MatDense<FPP,DEVICE>& U, const MatDense<FPP,DEVICE>& Lap, const ParamsFGFT<FPP,DEVICE,FPP2>& params, BlasHandle<DEVICE> cublasHandle, SpBlasHandle<DEVICE> cusparseHandle): HierarchicalFact<FPP, DEVICE, FPP2>(U, params, cublasHandle, cusparseHandle)
		{
			if ((U.getNbRow() != params.m_nbRow) |  (U.getNbCol() != params.m_nbCol))
				handleError(m_className,"constructor : params and Fourier matrix U haven't compatible size");
			if((Lap.getNbRow() != params.m_nbRow) |  (Lap.getNbCol() != params.m_nbCol))
				handleError(m_className,"constructor : params and Laplacian matrix Lap haven't compatible size");
			delete this->palm_global;
			this->palm_global = new Palm4MSAFGFT<FPP,DEVICE,FPP2>(Lap, params, cublasHandle, true);
		}


		const MatSparse<FPP, DEVICE>& get_D();

		void get_D(FPP* out_diag) const;

		void next_step()
		{
			this->palm_2.m_lambda = FPP2(1.);
			Faust::HierarchicalFact<FPP, DEVICE, FPP2>::next_step();
		}


	};
}

#include "faust_HierarchicalFactFGFT.hpp"

#endif
