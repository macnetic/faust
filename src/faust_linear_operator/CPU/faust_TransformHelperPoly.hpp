#include <cstdlib>
#include "faust_linear_algebra.h"
namespace Faust
{

	template<typename FPP>
	Vect<FPP,Cpu> TransformHelperPoly<FPP>::multiply(const Vect<FPP,Cpu> x, const bool transpose/*=false*/, const bool conjugate/*=false*/)
	{
		int d = L.getNbRow();
		Vect<FPP,Cpu> v1(d), v2, tmp;
		FPP *v0_buf = nullptr;
		int K = this->size();
		// seeing how basisChebyshev method is instantiating this class
		// K can't be strictly lower than two
		memcpy(v1.getData(), x.getData(), sizeof(FPP)*d);
		v2 = L*x;
		for(int i=2;i<=K;i++)
		{
			v0_buf = (FPP*) realloc(v0_buf, sizeof(FPP)*d*(i-1));
			memcpy(v0_buf+d*(i-2), v1.getData(), sizeof(FPP)*d);
			tmp = v1;
			v1 = v2;
//			gemv(L, v2, tmp, FPP(1), FPP(-1));
			v2 = twoL*v2;
			v2 -= tmp;
		}
		Vect<FPP, Cpu> v(d*(K+1));
		// K >= 2 so v0_buf is not nullptr
		memcpy(v.getData(), v0_buf, sizeof(FPP)*d*(K-1));
		memcpy(v.getData()+d*(K-1), v1.getData(), sizeof(FPP)*d);
		memcpy(v.getData()+d*K, v2.getData(), sizeof(FPP)*d);
		free(v0_buf);
		return v;
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K)
		{
			// assuming L is symmetric
			TransformHelper<FPP, Cpu>* basis = nullptr;
			MatSparse<FPP,Cpu> *T1, *T2;
			auto d = L->getNbRow();
			MatSparse<FPP,Cpu> Id, twoL, minus_Id, rR, R, zero;
			std::vector<MatGeneric<FPP,Cpu>*> facts(K);
			if(K == 0)
				return TransformHelper<FPP,Cpu>::eyeFaust(d,d);
			facts.resize(K); //normally it'd be K+1 but T0 is ignored
			// build the chebyshev polynomials by factor
			// K > 1 ignore the first one (for T0 0-degree), which is the identity
			// Identity
			Id.resize(d, d, d);
			Id.setEyes();
			// T1
			T1 = new MatSparse<FPP,Cpu>();
			T1->vstack(Id, *L);
			facts[K-1] = T1;
			if(K == 1)
				return new TransformHelper<FPP,Cpu>(facts, (FPP)1.0,
						/* optimizedCopy */ false,
						/* cloning_fact */ false,
						/* internal_call */ true);
			// rR block
			twoL = *L;
			twoL *= FPP(2);
			minus_Id.resize(d,d,d);
			minus_Id.setEyes();
			minus_Id *= FPP(-1);
			rR.hstack(minus_Id, twoL);
			// T2
			Id.resize(2*d, 2*d, 2*d);
			Id.setEyes();
			T2 = new MatSparse<FPP,Cpu>();
			T2->vstack(Id, rR);
			facts[K-2] = T2;
			// T3 to TK
			for(int i=3;i<K+1;i++)
			{
				Id.resize(i*d, i*d, i*d);
				Id.setEyes();
				zero.resize(0, d, (i-2)*d);
				R.hstack(zero, rR);
				auto Ti = new MatSparse<FPP,Cpu>();
				Ti->vstack(Id, R);
				facts[K-i] = Ti;
			}
			auto basisP = new TransformHelperPoly<FPP>(facts, (FPP)1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
			basisP->L = *L;
			basisP->twoL = *L;
			basisP->twoL *= 2;
			return basisP;
		}
}
