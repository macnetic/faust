#include <cstdlib>
#include "faust_linear_algebra.h"
namespace Faust
{

	template<typename FPP>
	Vect<FPP,Cpu> TransformHelperPoly<FPP>::multiply(const Vect<FPP,Cpu> x, const bool transpose/*=false*/, const bool conjugate/*=false*/)
	{

//		std::cout << "TransformHelperPoly<FPP>::multiply(Vect)" << std::endl;
		/**
		 * Recurrence relation (k=1 to K):
		 * v_{0,k} := [v_{0,k-1} ; v_{1, k-1} ] // concatenation
		 * v_{1,k} := v_{2,k-1}
		 * v_{2,k} := 2Lv_{2,k-1} - v_{1,k-1}
		 *
		 * First terms:
		 * v_{0,1} = []
		 * v_{1,1} = x
		 * v_{2,1} = Lx
		 *
		 * Fx = [v_{0, K} ; v_{1, K} ; v_{2, K}]
		 */
		int d = L.getNbRow();
		Vect<FPP,Cpu> v1(d), v2, tmp;
		int K = this->size();
		Vect<FPP, Cpu> v0(d*(K+1));
		// seeing how basisChebyshev method is instantiating this class
		// K can't be strictly lower than two
		memcpy(v1.getData(), x.getData(), sizeof(FPP)*d);
		v2 = L*x;
		for(int i=2;i<=K;i++)
		{
			memcpy(v0.getData()+d*(i-2), v1.getData(), sizeof(FPP)*d);
			tmp = v1;
			v1 = v2;
			//			gemv(L, v2, tmp, FPP(1), FPP(-1));
			v2 = twoL*v2;
			v2 -= tmp;
		}
		memcpy(v0.getData()+d*(K-1), v1.getData(), sizeof(FPP)*d);
		memcpy(v0.getData()+d*K, v2.getData(), sizeof(FPP)*d);
		return v0;
	}

// Slower method but kept commented here until further tests
//	template<typename FPP>
//			MatDense<FPP, Cpu> TransformHelperPoly<FPP>::multiply(MatDense<FPP,Cpu> X, const bool transpose/*=false*/, const bool conjugate/*=false*/)
//			{
////				std::cout << "TransformHelperPoly<FPP>::multiply(MatDense)" << std::endl;
//				int d = L.getNbRow();
//				int n = X.getNbCol();
//				int K = this->size();
//				MatDense<FPP, Cpu> V0(d*(K+1), n);
//				MatDense<FPP,Cpu> V1(d, n), V2 = X, tmp;
//				memcpy(V1.getData(), X.getData(), sizeof(FPP)*d*n);
//				L.multiply(V2, 'N');
//				for(int i=2;i<=K;i++)
//				{
//					memcpy(V0.getData()+d*n*(i-2), V1.getData(), sizeof(FPP)*d*n);
//					tmp = V1;
//					V1 = V2;
//					twoL.multiply(V2, 'N');
//					V2 -= tmp;
////					gemm(twoL, V2, tmp, (FPP)1.0, (FPP)-1.0, 'N', 'N');
////					V2  = tmp;
//				}
//				memcpy(V0.getData()+d*n*(K-1), V1.getData(), sizeof(FPP)*d*n);
//				memcpy(V0.getData()+d*n*K, V2.getData(), sizeof(FPP)*d*n);
//				// put bocks in proper order
//				MatDense<FPP,Cpu> V0_ord(d*(K+1), n);
//				#pragma omp parallel for
//				for(int i=0;i<(K+1);i++)
//				{
//					auto i_offset = i*n*d;
//					for(int j=0;j<n;j++)
//						memcpy(V0_ord.getData()+(K+1)*d*j+d*i, V0.getData()+j*d+i_offset, sizeof(FPP)*d);
//				}
//				return V0_ord;
//			}

	template<typename FPP>
			MatDense<FPP, Cpu> TransformHelperPoly<FPP>::multiply(MatDense<FPP,Cpu> X, const bool transpose/*=false*/, const bool conjugate/*=false*/)
			{
				//				std::cout << "TransformHelperPoly<FPP>::multiply(MatDense)" << std::endl;
				int d = L.getNbRow();
				int K = this->size();
				int n = X.getNbCol();
				MatDense<FPP,Cpu> V0_ord(d*(K+1), n);
				auto scale = (K+1)*d;
				#pragma omp parallel for
				for(int i=0;i<n;i++)
				{
					Vect<FPP,Cpu> x(d, X.getData()+i*d);
					auto y = multiply(x);
					memcpy(V0_ord.getData()+scale*i, y.getData(), sizeof(FPP)*scale);
				}
				return V0_ord;
			}

	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperPoly<FPP>::poly(MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs)
		{
			coeffs.multiplyLeft(basisX);
			return coeffs;
		}

	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperPoly<FPP>::poly(int n, MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs)
		{
			int d = L.getNbRow();
			MatDense<FPP,Cpu> Y(d, n);
			this->poly(d, n, basisX.getData(), coeffs.getData(), Y.getData());
			return Y;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::poly(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			int K = this->size();
			Faust::poly(d, K, n, basisX, coeffs, out);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::polyFaust(const FPP* coeffs)
		{
//            Id = sp.eye(L.shape[1], format="csr")
//				            scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
//									                                format="csr")
			MatSparse<FPP,Cpu> Id, Id1;
			MatSparse<FPP,Cpu> coeffDiags;
			auto d = L.getNbRow();
			int K = this->size();
			Id.resize(d, d, d);
			Id.setEyes();
			auto Id0 = Id;
			Id0 *= coeffs[0];
			// K >= 2
			for(int i=1;i < K+1; i++)
			{
				Id1 = Id;
				Id1 *= coeffs[i];
				MatSparse<FPP,Cpu> coeffDiags;
				coeffDiags.hstack(Id0, Id1);
				Id0 = coeffDiags;
			}
			std::vector<MatGeneric<FPP,Cpu>*> facts(K+1);
			facts.resize(K+1);
			facts[0] = new MatSparse<FPP,Cpu>(Id0);
			for(int i=1;i<=K;i++)
			{
				facts[i] = this->get_gen_fact_nonconst(i-1);
			}
			return new TransformHelper<FPP,Cpu>(facts, (FPP) 1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K)
		{
			// assuming L is symmetric
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
			#pragma omp parallel for private(Id, zero, R)
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

	template<typename FPP>
		void poly(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			auto K_plus_1 = K+1;
			auto d_K_plus_1 = d*K_plus_1;
			#pragma omp parallel for
			for(int i=0;i<n;i++)
			{
				Vect<FPP,Cpu> _coeffs(K_plus_1, coeffs);
				MatDense<FPP,Cpu> basisXi(d, K_plus_1);
				memcpy(basisXi.getData(), basisX+d_K_plus_1*i, sizeof(FPP)*d_K_plus_1);
				_coeffs.multiplyLeft(basisXi);
				memcpy(out+i*d, _coeffs.getData(), sizeof(FPP)*d);
			}
		}
}
