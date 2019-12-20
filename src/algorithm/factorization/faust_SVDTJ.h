#ifndef __FAUST_SVDTJ__
#define __FAUST_SVDTJ__
#include "faust_MatDense.h"
#include "faust_MatSparse.h"

#include "faust_constant.h"
namespace Faust
{
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj(MatDense<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj(MatSparse<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj_cplx(MatDense<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj_cplx(MatSparse<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);

template<typename FPP, Device DEVICE, typename FPP2 = float>
void svdtj_core_gen(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, double tol, unsigned int verbosity, bool relErr, int order /* ignored anyway, kept here just in case */, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_);

};

#include "faust_SVDTJ.hpp"
#endif
