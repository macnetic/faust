#include "faust_MatDense.h"
#include "faust_MatSparse.h"

namespace Faust
{
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj(MatDense<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj(MatSparse<FPP, DEVICE> & M, int J, int t, double tol, unsigned int verbosity, bool relErr, int order,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S);
	template<typename FPP, Device DEVICE, typename FPP2 = float>
		void svdtj_core(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_);

};

#include "faust_SVDTJ.hpp"
