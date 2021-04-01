#include <cstdlib>
#include "faust_linear_algebra.h"
#include "Eigen/Core"

namespace Faust
{
	template<typename FPP>
		faust_unsigned_int TransformHelperPoly<FPP>::getNbRow() const
		{
			int d = L->getNbRow();
			int nrows, ncols;
			nrows = this->size()*d; // (K+1)*d
			ncols = d;
			if(this->is_sliced){
				faust_unsigned_int id = this->is_transposed?1:0;
				return	this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return this->is_transposed?ncols:nrows;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelperPoly<FPP>::getNbCol() const
		{
			int d = L->getNbRow();
			int nrows, ncols;
			nrows = this->size()*d; // (K+1)*d
			ncols = d;
			if(this->is_sliced)
			{
				faust_unsigned_int id = this->is_transposed?0:1;
				return this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return this->is_transposed?nrows:ncols;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelperPoly<FPP>::multiply(const Vect<FPP,Cpu> &x, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			return std::move(this->multiply(x.getData(), transpose, conjugate));
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelperPoly<FPP>::multiply(const FPP* x, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			int d = L->getNbRow();
			int K = this->size()-1;
			Vect<FPP, Cpu> v0(d*(K+1));
			multiply(x, v0.getData(), transpose, conjugate);
			return std::move(v0);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply(const FPP* x, FPP* y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
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
			int d = L->getNbRow();
			int K = this->size()-1;

			memcpy(y, x, sizeof(FPP)*d);
			if(K == 0)
				return;
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> x_vec(const_cast<FPP*>(x), d);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>v2(const_cast<FPP*>(y+d), d);
			v2 = L->mat*x_vec;
			if(K == 1) // not necessary but clearer
				return;
			for(int i=3;i<=K+1;i++)
			{
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>new_v2_(const_cast<FPP*>(y+d*(i-1)), d);
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>v2_(const_cast<FPP*>(y+d*(i-2)), d);
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>v1_(const_cast<FPP*>(y+d*(i-3)), d);
				new_v2_ = L->mat*v2_*2-v1_;
			}
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
		MatDense<FPP, Cpu> TransformHelperPoly<FPP>::multiply(const MatDense<FPP,Cpu> &X, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			//				std::cout << "TransformHelperPoly<FPP>::multiply(MatDense)" << std::endl;
			int d = L->getNbRow();
			int K = this->size()-1;
			int n = X.getNbCol();
			MatDense<FPP,Cpu> Y(d*(K+1), n);
			multiply(X.getData(), n, Y.getData(), transpose, conjugate);
			return Y;
		}

template<typename FPP>
		void TransformHelperPoly<FPP>::multiply(const FPP* X, int n, FPP* Y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			int d = L->getNbRow();
			int K = this->size()-1;
			auto scale = (K+1)*d;
#pragma omp parallel for
			for(int i=0;i<n;i++)
			{
				//				Vect<FPP,Cpu> x(d, X.getData()+i*d);
				//				auto y = multiply(x);
				//				memcpy(V0_ord.getData()+scale*i, y.getData(), sizeof(FPP)*scale);
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> x_vec(const_cast<FPP*>(X+i*d), d);
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> y_vec(const_cast<FPP*>(Y+i*scale), scale);
				multiply(x_vec.data(), y_vec.data(), transpose, conjugate);
			}
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
			int d = L->getNbRow();
			MatDense<FPP,Cpu> Y(d, n);
			this->poly(d, n, basisX.getData(), coeffs.getData(), Y.getData());
			return Y;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::poly(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			int K = this->size()-1;
			Faust::poly(d, K, n, basisX, coeffs, out);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::next()
		{
			int K = this->size()-1;
			return next(K+1);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::next(int K)
		{
			auto old_K = this->size()-1;
			auto d = L->getNbRow();
			MatSparse<FPP,Cpu> Id, zero, R;
			std::vector<MatGeneric<FPP,Cpu>*> facts(K+1);
			for(int i=0;i<old_K+1;i++)
				facts[K-i] = this->get_gen_fact_nonconst(old_K-i);
#pragma omp parallel for private(Id, zero, R)
			for(int i=old_K+1;i <= K; i++)
			{
				//TODO: refactor with basisChebyshev
				Id.resize(i*d, i*d, i*d);
				Id.setEyes();
				zero.resize(0, d, (i-2)*d);
				R.hstack(zero, *rR);
				auto Ti = new MatSparse<FPP,Cpu>();
				Ti->vstack(Id, R);
				facts[K-i] = Ti;
			}
			auto basisP = new TransformHelperPoly<FPP>(facts, (FPP)1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
			basisP->L = L;
			ref_man.acquire(L);
			basisP->rR = rR;
			ref_man.acquire(rR);
			return basisP;
		}


	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::polyFaust(const FPP* coeffs)
		{
			//            Id = sp.eye(L.shape[1], format="csr")
			//				            scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
			//									                                format="csr")
			std::vector<MatGeneric<FPP,Cpu>*> facts(this->size()+1);
			MatSparse<FPP,Cpu> Id, Id1;
			MatSparse<FPP,Cpu> coeffDiags, tmp;
			auto d = L->getNbRow();
			Id.resize(d, d, d);
			Id.setEyes();
			coeffDiags = Id;
//			copy_sp_mat(Id, coeffDiags);
			coeffDiags *= coeffs[0];
			for(int i=1;i < this->size(); i++)
			{
				Id1 = Id;
//				copy_sp_mat(Id, Id1);
				Id1 *= coeffs[i];
				tmp.hstack(coeffDiags, Id1); //TODO: hstack/vstack to concatenate "this" to argument matrix
//				coeffDiags = tmp; // crash //TODO: fix
				copy_sp_mat(tmp, coeffDiags);
			}
			coeffDiags.set_id(false);
			auto fac0 = new MatSparse<FPP,Cpu>(); // ctor copy crash //TODO: fix
			copy_sp_mat(coeffDiags, *fac0);
			facts[0] = fac0;
			for(int i=1;i <= this->size();i++)
				facts[i] = this->get_gen_fact_nonconst(i-1);
			return new TransformHelper<FPP,Cpu>(facts, (FPP) 1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
		}

	template<typename FPP>
		Faust::RefManager TransformHelperPoly<FPP>::ref_man([](void *fact)
				{
#ifdef DEBUG
				std::cout << "Faust::TransformHelperPoly delete" << std::endl;
#endif
				delete static_cast<MatGeneric<FPP,Cpu>*>(fact);
				});

	template<typename FPP>
	TransformHelperPoly<FPP>::~TransformHelperPoly()
	{
#ifdef FAUST_VERBOSE
		std::cout << "~TransformHelperPoly()" << std::endl;
#endif
		ref_man.release(L);
		ref_man.release(rR);
	}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshevT0(MatSparse<FPP,Cpu>* T0/*=nullptr*/)
		{
			int K = this->size()-1;
			if(T0 == nullptr)
			{
				// Identity
				auto T0_ = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(K));
				auto d = this->L->getNbRow();
				T0_->resize(d, d, d);
				T0_->setEyes();
			}
			else
				this->update(*T0, K);
			is_fact_created[K] = true;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshevT1()
		{
			int K = this->size()-1;
			MatSparse<FPP,Cpu> Id;
			int id = K-1;
			if(! is_fact_created[id])
			{
				auto d = this->L->getNbRow();
				Id.resize(d, d, d);
				Id.setEyes();
				auto T1 = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(id));
				T1->vstack(Id, *L);
				is_fact_created[id] = true;
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshevT2()
		{
			int K = this->size()-1;
			int id = K-2;
			if(! is_fact_created[id])
			{
				auto d = this->L->getNbRow();
				MatSparse<FPP,Cpu> Id;
				Id.resize(2*d, 2*d, 2*d);
				Id.setEyes();
				auto T2 = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(id));
				T2->vstack(Id, *rR);
				is_fact_created[id] = true;
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshevTi(int i)
		{
			switch(i)
			{
				case 0:
					return basisChebyshevT0();
				case 1:
					return basisChebyshevT1();
				case 2:
					return basisChebyshevT2();
				default:
					int K = this->size()-1;
					int fid = K-i;
					if(! is_fact_created[fid])
					{
						MatSparse<FPP,Cpu> R, zero;
						MatSparse<FPP,Cpu> Id;
						auto d = this->L->getNbRow();
						auto id = i*d;
						Id.resize(id, id, id);
						Id.setEyes();
						zero.resize(d, (i-2)*d);
						R.hstack(zero, *rR);
						auto Ti = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(fid));
						Ti->vstack(Id, R);
						is_fact_created[fid] = true;
					}
					break;
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_all()
		{
			int K = this->size()-1;
			this->basisChebyshevT0();

			if (K > 0)
			{
				this->basisChebyshevT1();
				if (K > 1)
				{
					this->basisChebyshevT2();
					#pragma omp parallel for
					for(int i=3;i<K+1;i++)
						this->basisChebyshevTi(i);
				}
			}

		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::create_rR(const MatSparse<FPP,Cpu>* L)
		{
			MatSparse<FPP,Cpu> twoL, minus_Id;
			auto d = this->L->getNbRow();
			// rR block
			twoL = *L;
			twoL *= FPP(2);
			minus_Id.resize(d,d,d);
			minus_Id.setEyes();
			minus_Id *= FPP(-1);
			rR = new MatSparse<FPP,Cpu>();
			rR->hstack(minus_Id, twoL);
			this->rR = rR;
			Faust::TransformHelperPoly<FPP>::ref_man.acquire(rR);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0/*=nullptr*/, bool lazy_instantiation/*=true*/)
		{

			// assuming L is symmetric
			MatSparse<FPP,Cpu> *T1, *T2;
			MatSparse<FPP,Cpu> Id, twoL, minus_Id, *rR, R, zero;
			auto d = L->getNbRow();

			std::vector<MatGeneric<FPP,Cpu>*> facts(K+1);
			for(int i=0;i<K+1;i++)
			{
				facts[i] = new MatSparse<FPP,Cpu>(); // empty factors (lazy instantiation)
			}
			auto basisP = new TransformHelperPoly<FPP>(facts, (FPP)1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
			basisP->is_fact_created.assign(K+1, false);
			basisP->L = new MatSparse<FPP,Cpu>(*L);
			Faust::TransformHelperPoly<FPP>::ref_man.acquire(basisP->L);

			basisP->create_rR(L);

			// build the chebyshev polynomials by factor
			if(T0 != nullptr)
				// T0 is arbitrary and need to be stored
				// best choice is to instantiate the factor directly
				basisP->basisChebyshevT0(T0);

			if(! lazy_instantiation)
				basisP->basisChebyshev_all();
			return basisP;
		}

	template<typename FPP>
		void poly(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			auto K_plus_1 = K+1;
			auto d_K_plus_1 = d*K_plus_1;
			//			#pragma omp parallel for // most likely counterproductive to use OpenMP here, Eigen is already multithreading scalar product in the loop
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> vec_out(out, d, n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> vec_coeffs(const_cast<FPP*>(coeffs), K_plus_1);
			for(int i=0;i<n;i++)
			{
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> mat_basisX(const_cast<FPP*>(basisX+d_K_plus_1*i), d, K_plus_1); //constcast is not dangerous, no modification is done later
				vec_out.block(0,i,d,1) = mat_basisX*vec_coeffs;
			}
		}
}
