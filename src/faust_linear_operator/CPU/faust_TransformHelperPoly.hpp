#include <cstdlib>
#include "faust_linear_algebra.h"
#include "Eigen/Core"
#define FUI faust_unsigned_int

namespace Faust
{
	template<typename FPP>
		TransformHelperPoly<FPP>::TransformHelperPoly(uint K,
				MatSparse<FPP, Cpu>* L,
				MatSparse<FPP, Cpu>* rR /* = nullptr*/,
				MatSparse<FPP,Cpu> *T0 /*= nullptr*/,
				BasisLaziness laziness /*= INSTANTIATE_COMPUTE_AND_FREE */,
				bool on_gpu /*= false*/) : TransformHelper<FPP,Cpu>()
	{
		// assuming L is symmetric
		this->L = L;
		uint d = L->getNbRow();
		ref_man.acquire(L);
		this->rR = rR;
		if(rR == nullptr)
			this->create_rR(L);
		ref_man.acquire(rR);
		this->laziness = laziness;

		this->is_fact_created.assign(K+1, laziness == NOT_LAZY);
		// build the chebyshev polynomials by factor
		for(int i=0;i<K+1;i++)
		{
			// empty factors (in case of lazy instantiation)
			this->push_back(new MatSparse<FPP, Cpu>()/*d*(K-i+1), d*(K-i)*/, /* OptimizedCopy */ false, /* copying */ false, /*transpose*/ false, /*conjugate*/ false, /* verify_dims_agree */ false); // the dimensions are fake
		}

		if(T0 != nullptr)
			// T0 is arbitrary and need to be stored
			// best choice is to instantiate the factor directly
			this->basisChebyshevT0(T0);
			//T0_is_arbitrary == true
			//ref_man acquired T0

		if(laziness == NOT_LAZY)
			this->basisChebyshev_all();

		this->mul_and_combi_lin_on_gpu = on_gpu;
	}

	template<typename FPP>
		TransformHelperPoly<FPP>::TransformHelperPoly(uint K, const TransformHelperPoly<FPP>& src) : TransformHelper<FPP, Cpu>()
	{
		if(K+1 < src.size())
			throw std::runtime_error("The src TransformHelperPoly size can't be greater than K+1.");
		uint d = src.L->getNbRow();
		int i;
		int num_new_polys = (K+1)-src.size();
		for(i=0;i<num_new_polys;i++)
		{
			this->push_back(new MatSparse<FPP, Cpu>()/*d*(K-i+1), d*(K-i)*/, /* OptimizedCopy */ false, /* copying */ false, /*transpose*/ false, /*conjugate*/ false, /* verify_dims_agree */ false);
		}
		for(;i<src.size()+num_new_polys;i++)
			this->push_back(src.get_gen_fact_nonconst(i-num_new_polys), /* OptimizedCopy */ false, /* copying */ false, /*transpose*/ false, /*conjugate*/ false, /* verify_dims_agree */ false);
		// copy all attributes
		laziness = src.laziness;
		this->is_fact_created.assign(this->size(), laziness == NOT_LAZY);
		for(i=num_new_polys; i < this->size();i++)
		{
			is_fact_created[i] = is_fact_created[i-num_new_polys];
		}
		L = src.L;
		ref_man.acquire(L);
		rR = src.rR;
		ref_man.acquire(rR);
		T0_is_arbitrary = src.T0_is_arbitrary;
		if(laziness == NOT_LAZY)
			this->basisChebyshev_facti2j(0, K);
	}


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
			uint K = this->size()-1;
			Vect<FPP, Cpu> v0(d*(K+1));
			multiply(x, v0.getData(), transpose, conjugate);
			return std::move(v0);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply_cpu(const FPP* x, FPP* y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
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
			uint K = this->size()-1;

			memcpy(y, x, sizeof(FPP)*d);
			if(K == 0)
				return;
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> x_vec(const_cast<FPP*>(x), d);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> v2(const_cast<FPP*>(y+d), d);
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

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply_gpu(const FPP* x, FPP* y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
#ifdef USE_GPU_MOD
			int d = L->getNbRow();
			uint K = this->size()-1;
			Vect<FPP, GPU2> gpu_v1(d, x);
			Vect<FPP, GPU2> gpu_v2(gpu_v1);
			Vect<FPP, GPU2> gpu_new_v2(d);
			const MatSparse<FPP, GPU2> gpu_L(*this->L);
			MatSparse<FPP, GPU2> gpu_twoL(gpu_L);
			gpu_twoL *= 2;
			memcpy(y, x, sizeof(FPP)*d);
			if(K == 0)
				return;
			//			gpu_v2 == x
			gpu_v2.multiplyLeft(gpu_L);
			gpu_v2.tocpu(y+d); //			v2 = L->mat*x_vec;
			if(K == 1) // not necessary but clearer
				return;
			for(int i=3;i<=K+1;i++)
			{
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>new_v2_cpu(const_cast<FPP*>(y+d*(i-1)), d);
				gpu_new_v2 = gpu_v2;
				gpu_new_v2.multiplyLeft(const_cast<const MatSparse<FPP,GPU2>&>(gpu_twoL));
				gpu_new_v2 -= gpu_v1; // new_v2_ = L->mat*v2_*2-v1_;
				// prepare next it
				gpu_v1 = gpu_v2;
				gpu_v2 = gpu_new_v2;
				gpu_new_v2.tocpu(new_v2_cpu.data());
			}
#else
			throw std::runtime_error("USE_GPU_MOD option must be enabled at compiling time to use this function (TransformHelperPoly<FPP>::multiply_gpu).");
#endif
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply(const FPP* x, FPP* y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			if(this->mul_and_combi_lin_on_gpu)
			{
				multiply_gpu(x, y, transpose, conjugate);
			}
			else
			{
				multiply_cpu(x, y, transpose, conjugate);
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly(const FPP* x, FPP* y, const FPP* coeffs)
		{
			if(this->mul_and_combi_lin_on_gpu)
				this->multiplyPoly_gpu(x, y, coeffs);
			else
				this->multiplyPoly_cpu(x, y, coeffs);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly_cpu(const FPP* x, FPP* y, const FPP* coeffs)
		{
			multiplyPoly_cpu(x, 1, y, coeffs);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly_gpu(const FPP* x, FPP* y, const FPP* coeffs)
		{
			multiplyPoly_gpu(x, 1, y, coeffs);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly(const FPP* X, int n, FPP* Y, const FPP* coeffs)
		{
			if(this->mul_and_combi_lin_on_gpu)
				multiplyPoly_gpu(X, n, Y, coeffs);
			else
				multiplyPoly_cpu(X, n, Y, coeffs);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly_cpu(const FPP* X, int n, FPP* Y, const FPP* coeffs)
		{
			int d = L->getNbRow();
			uint K = this->size()-1;
			memcpy(Y, X, sizeof(FPP)*d*n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> y_(const_cast<FPP*>(Y), d, n);
			y_ *= coeffs[0];
			if(K == 0)
				return;
			FPP* buf = new FPP[3*d*n]; // working buffer for v0, v1, v2
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> x_vec(const_cast<FPP*>(X), d, n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> v0(const_cast<FPP*>(buf), d, n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> v1(const_cast<FPP*>(buf+d*n), d, n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> v2(const_cast<FPP*>(buf+2*d*n), d, n);
			memcpy(v0.data(), X, sizeof(FPP)*d*n);
			v1 = L->mat * x_vec;
			y_ += v1*coeffs[1];
			if(K == 1)
				return;
			for(int i=2;i<=K;i++)
			{
				v2 = L->mat * v1 * 2 - v0;
				y_ += v2*coeffs[i];
				v0 = v1;
				v1 = v2;
//				memcpy(v0.data(), v1.data(), sizeof(FPP)*d*n);
//				memcpy(v1.data(), v2.data(), sizeof(FPP)*d*n);
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiplyPoly_gpu(const FPP* X, int n, FPP* Y, const FPP* coeffs)
		{
#ifdef USE_GPU_MOD
			int d = L->getNbRow();
			uint K = this->size()-1;
			MatDense<FPP, GPU2> gpu_V1((FUI) d, (FUI) n, X);
			MatDense<FPP, GPU2> gpu_V2(gpu_V1);
			MatDense<FPP, GPU2> gpu_new_V2((FUI) d, (FUI) n);
			MatDense<FPP, GPU2> gpu_Y((FUI) d, (FUI) n, X);
			const MatSparse<FPP, GPU2> gpu_L(*this->L);
			MatSparse<FPP, GPU2> gpu_twoL(gpu_L);
			gpu_twoL *= 2;
			gpu_Y.scalarMultiply(coeffs[0]); // x*coeffs[0]
			if(K == 0)
				return;
			//			gpu_V2 == x
			gpu_V2.multiplyLeft(gpu_L);
			gpu_new_V2 = gpu_V2;
			gpu_new_V2.scalarMultiply(coeffs[1]); //	coeffs[1]*(L->mat*x);
			gpu_Y.add(gpu_new_V2); // gpu_Y = x*coeffs[0]+coeffs[1]*(L->mat*x)
			if(K == 1) // not necessary but clearer
				return;
			for(int i=3;i<=K+1;i++)
			{
				gpu_new_V2 = gpu_V2;
				gpu_new_V2.multiplyLeft(const_cast<const MatSparse<FPP,GPU2>&>(gpu_twoL));
				gpu_new_V2 -= gpu_V1; // new_v2_ = L->mat*v2_*2-v1_;
				// prepare next it
				gpu_V1 = gpu_V2;
				gpu_V2 = gpu_new_V2;
				gpu_new_V2.scalarMultiply(coeffs[i-1]);
				gpu_Y.add(gpu_new_V2);
			}
			gpu_Y.tocpu(Y, nullptr); // warning: without explicitely passing a nullptr for the stream it provokes a cudaErrorUnknown (TODO: verify tocpu overloads)
#else
			throw std::runtime_error("USE_GPU_MOD option must be enabled at compiling time to use this function (TransformHelperPoly<FPP>::multiplyPoly_gpu).");
#endif
		}

	// Slower method but kept commented here until further tests
	//	template<typename FPP>
	//			MatDense<FPP, Cpu> TransformHelperPoly<FPP>::multiply(MatDense<FPP,Cpu> X, const bool transpose/*=false*/, const bool conjugate/*=false*/)
	//			{
	////				std::cout << "TransformHelperPoly<FPP>::multiply(MatDense)" << std::endl;
	//				int d = L.getNbRow();
	//				int n = X.getNbCol();
	//				uint K = this->size();
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
			uint K = this->size()-1;
			int n = X.getNbCol();
			MatDense<FPP,Cpu> Y(d*(K+1), n);
			multiply(X.getData(), n, Y.getData(), transpose, conjugate);
			return Y;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply(const FPP* X, int n, FPP* Y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			if(this->mul_and_combi_lin_on_gpu)
			{
				multiply_gpu(X, n, Y, transpose, conjugate);
			}
			else
			{
				multiply_cpu(X, n, Y, transpose, conjugate);
			}
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::multiply_cpu(const FPP* X, int n, FPP* Y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			int d = L->getNbRow();
			uint K = this->size()-1;
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
		void TransformHelperPoly<FPP>::multiply_gpu(const FPP* X, int n, FPP* Y, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
#ifdef USE_GPU_MOD
			int d = L->getNbRow();
			uint K = this->size()-1;
			MatDense<FPP, GPU2> gpu_V1((FUI) d, (FUI) n, X);
			MatDense<FPP, GPU2> gpu_V2(gpu_V1);
			MatDense<FPP, GPU2> gpu_new_V2((FUI) d, (FUI) n);
			MatDense<FPP, Cpu> tmp_cpu_V2((FUI)d, (FUI)n);
			const MatSparse<FPP, GPU2> gpu_L(*this->L);
			MatSparse<FPP, GPU2> gpu_twoL(gpu_L);
			gpu_twoL *= 2;
			auto block_to_cpu = [&Y, &d, &n, &K, &tmp_cpu_V2](int i, const FPP* X_)
			{
				#pragma omp parallel for
				for(int j=0;j<n;j++)
				{
					memcpy(Y+(K+1)*d*j+d*i, X_+j*d, sizeof(FPP)*d);
				}
			};
			block_to_cpu(0, X);
			if(K == 0)
				return;
			//			gpu_V2 == X
			gpu_V2.multiplyLeft(gpu_L);
			gpu_V2.tocpu(tmp_cpu_V2); //			v2 = L->mat*x_vec;
			block_to_cpu(1, tmp_cpu_V2.getData());
			if(K == 1) // not necessary but clearer
				return;
			for(int i=3;i<=K+1;i++)
			{
				gpu_new_V2 = gpu_V2;
				gpu_new_V2.multiplyLeft(const_cast<const MatSparse<FPP,GPU2>&>(gpu_twoL));
				gpu_new_V2 -= gpu_V1; // new_v2_ = L->mat*v2_*2-v1_;
				// prepare next it
				gpu_V1 = gpu_V2;
				gpu_V2 = gpu_new_V2;
				gpu_new_V2.tocpu(tmp_cpu_V2);
				block_to_cpu(i-1, tmp_cpu_V2.getData());
			}
#else
			throw std::runtime_error("USE_GPU_MOD option must be enabled at compiling time to use this function (TransformHelperPoly<FPP>::multiply_gpu).");
#endif
		}

	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperPoly<FPP>::multiply(const MatSparse<FPP,Cpu> &A, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			MatDense<FPP, Cpu> M(this->getNbRow(), A.getNbCol());
			M.setZeros();
			vector<int> col_ids;
			for(int i=0;i<A.getNonZeros();i++)
			{
				auto id = A.getColInd()[i];
				if(std::find(std::begin(col_ids), std::end(col_ids), id) == std::end(col_ids))
					col_ids.push_back(id);
			}
#ifndef _MSC_VER // MS VS doesn't want this OMP loop (increment error) // error C3017
			#pragma omp parallel for
#endif
			for(auto i=col_ids.begin(); i < col_ids.end();i++)
			{
				auto id = *i;
				auto vect = A.get_col(id);
				this->multiply(vect.getData(), M.getData()+M.getNbRow()*id, transpose, conjugate);
			}
			return M;
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
			uint K = this->size()-1;
			Faust::poly(d, K, n, basisX, coeffs, out, this->mul_and_combi_lin_on_gpu);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::next()
		{
			uint K = this->size()-1;
			return next(K+1);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::next(uint K)
		{
			auto basisP = new TransformHelperPoly<FPP>(K, *this);
			return basisP;
		}


	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::polyFaust(const FPP* coeffs)
		{
			//            Id = sp.eye(L.shape[1], format="csr")
			//				            scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
			//									                                format="csr")
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
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
			auto ret = new TransformHelper<FPP,Cpu>(facts, (FPP) 1.0,
					/* optimizedCopy */ false,
					/* cloning_fact */ false,
					/* internal_call */ true);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
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
			uint K = this->size()-1;
			if(T0 == nullptr)
			{
				// Identity
				auto T0_ = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(K));
				auto d = this->L->getNbRow();
				T0_->resize(d, d, d);
				T0_->setEyes();
				T0_is_arbitrary = false;
			}
			else
			{
				this->update(*T0, K);
				T0_is_arbitrary = true;
			}
			is_fact_created[K] = true;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshevT1()
		{
			uint K = this->size()-1;
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
			uint K = this->size()-1;
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
		void TransformHelperPoly<FPP>::basisChebyshevTi(uint i)
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
					uint K = this->size()-1;
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
			uint K = this->size()-1;
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
		TransformHelper<FPP,Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0/*=nullptr*/, bool on_gpu/*=false*/, BasisLaziness lazy_instantiation/*=INSTANTIATE_COMPUTE_AND_FREE*/)
		{
			// assuming L is symmetric
			auto basisP = new TransformHelperPoly<FPP>(
					K,
					new MatSparse<FPP, Cpu>(*L),
					nullptr /* rR initialized by the ctor */,
					T0,
					lazy_instantiation,
					on_gpu);

			return basisP;
		}

	template<typename FPP>
		void poly(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out, bool on_gpu/*=on_gpu*/)
		{
			if(on_gpu)
#ifdef USE_GPU_MOD
				poly_gpu(d, K, n, basisX, coeffs, out);
#else
			throw std::runtime_error("USE_GPU_MOD option must be enabled at compiling time to use this function (Faust::poly_gpu()).");
#endif
			else // cpu
				poly_cpu(d, K, n, basisX, coeffs, out);
		}

	template<typename FPP>
		void poly_cpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			auto K_plus_1 = K+1;
			auto d_K_plus_1 = d*K_plus_1;
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> vec_out(out, d, n);
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> vec_coeffs(const_cast<FPP*>(coeffs), K_plus_1);
			//			#pragma omp parallel for // most likely counterproductive to use OpenMP here, Eigen is already multithreading scalar product in the loop
			for(int i=0;i<n;i++)
			{
				Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> mat_basisX(const_cast<FPP*>(basisX+d_K_plus_1*i), d, K_plus_1); //constcast is not dangerous, no modification is done later
				vec_out.block(0,i,d,1) = mat_basisX*vec_coeffs;
			}
		}


	template<typename FPP>
		void poly_gpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out)
		{
			auto K_plus_1 = K+1;
			auto d_K_plus_1 = d*K_plus_1;
			Vect<FPP, GPU2> gpu_vec_coeffs(K_plus_1, coeffs);
			Vect<FPP, GPU2> gpu_dsize_vec(d);
			for(int i=0;i<n;i++)
			{
				MatDense<FPP,GPU2> gpu_mat_basis((FUI)d, (FUI)K_plus_1, basisX+d_K_plus_1*i);
//				gpu_mat_basis.multiply(gpu_vec_coeffs, gpu_dsize_vec);
//				gpu_dsize_vec.tocpu(out+d*i);
				poly_gpu(gpu_mat_basis, gpu_vec_coeffs, gpu_dsize_vec, out, d*i);
			}
		}

	template<typename FPP>
		void poly_gpu(const MatDense<FPP,GPU2>& gpu_mat_basis, const Vect<FPP, GPU2>& gpu_vec_coeffs, Vect<FPP, GPU2>& gpu_dsize_vec, FPP* out, int out_offset)
		{
			gpu_mat_basis.multiply(gpu_vec_coeffs, gpu_dsize_vec);
			gpu_dsize_vec.tocpu(out+out_offset);
		}

	template<typename FPP>
		void poly(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out, bool on_gpu/*=false*/)
		{
			if(on_gpu)
#ifdef USE_GPU_MOD
				poly_gpu(d, K, n, basisX, coeffs, out, n_out);
#else
			throw std::runtime_error("USE_GPU_MOD option must be enabled at compiling time to use this function (Faust::poly_gpu()).");
#endif
			else // cpu
				poly_cpu(d, K, n, basisX, coeffs, out, n_out);

		}

	template<typename FPP>
		void poly_cpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out)
		{
			for(int i=0;i<n_out;i++)
				poly_cpu(d, K, n, basisX, coeffs+i*(K+1), out[i]);
		}

	template<typename FPP>
		void poly_gpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out)
		{
			auto K_plus_1 = K+1;
			auto d_K_plus_1 = d*K_plus_1;
			Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>> vec_coeffs(const_cast<FPP*>(coeffs), K_plus_1);
			Vect<FPP, GPU2>* gpu_vec_coeffs[n_out];
			Vect<FPP, GPU2> gpu_dsize_vec(d);
			for(int j=0;j<n_out;j++)
				gpu_vec_coeffs[j] = new Vect<FPP, GPU2>(K_plus_1, coeffs+j*(K+1));
			for(int i=0;i<n;i++)
			{
				MatDense<FPP,GPU2> gpu_mat_basis((FUI)d, (FUI)K_plus_1, basisX+d_K_plus_1*i);
				for(int j=0;j<n_out;j++)
				{
					poly_gpu(gpu_mat_basis, *gpu_vec_coeffs[j], gpu_dsize_vec, out[j], d*i);
				}
			}
			for(int j=0;j<n_out;j++)
				delete gpu_vec_coeffs[j];
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_facti(uint i)
		{
			uint K = this->size()-1; // size() >= 1 (cf. ctor)
			int k_i = K-i;
			basisChebyshevTi(k_i);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_free_facti(uint i)
		{
			if(i >= this->size()) throw std::out_of_range("i is greater than size");
			uint K = this->size()-1;

			if(i != 0 || ! T0_is_arbitrary)
			{
				auto Ti = dynamic_cast<MatSparse<FPP,Cpu>*>(this->get_gen_fact_nonconst(i));
				Ti->resize(0,0,0);
				is_fact_created[i] = false;
			}
			// else don't remove a T0 which comes from the oustide (we don't know how to recreate it later)
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_facti2j(uint i, uint j)
		{
			if(i > j) throw std::runtime_error("i must be lower than j");
			if(j >= this->size()) throw std::out_of_range("j is greater than size");
			for(auto k=i;k<=j;k++)
				basisChebyshev_facti(k);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_free_facti2j(uint i, uint j)
		{
			if(i > j) throw std::runtime_error("i must be lower than j");
			for(auto k=i;k<=j;k++)
				basisChebyshev_free_facti(k);
		}


	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_fact_all()
		{
			basisChebyshev_facti2j(0, this->size()-1);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::basisChebyshev_free_fact_all()
		{

			basisChebyshev_free_facti2j(0, this->size()-1);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::get_fact(const faust_unsigned_int &id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose/*= false*/) const
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_facti(id);
			TransformHelper<FPP, Cpu>::get_fact(id, elts, num_rows, num_cols, transpose);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_facti(id);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::get_fact(const faust_unsigned_int id,
				int* rowptr,
				int* col_ids,
				FPP* elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose/*= false*/) const
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_facti(id);
			TransformHelper<FPP, Cpu>::get_fact(id, rowptr, col_ids, elts, nnz, num_rows, num_cols, transpose);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_facti(id);
		}

	template<typename FPP>
		uint TransformHelperPoly<FPP>::get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const
		{
			//dim == 0 to get num of cols, >0 otherwise
			faust_unsigned_int rid; //real id
			if(this->is_transposed) {
				rid = this->size()-id-1;
				dim = !dim;
			}
			else {
				rid = id;
				//dim = dim;
			}
			if(dim)
				return this->L->getNbCol()*(this->size()-rid-(rid==this->size()-1?0:1));
			else
				return this->L->getNbCol()*(this->size()-rid);
		}

	template<typename FPP>
		const MatGeneric<FPP,Cpu>* TransformHelperPoly<FPP>::get_gen_fact(const faust_unsigned_int id) const
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_facti(id);
			auto ret = TransformHelper<FPP, Cpu>::get_gen_fact(id);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_facti(id);
			return ret;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelperPoly<FPP>::get_fact_nnz(const faust_unsigned_int id) const
		{
			// compute the nnz based on the factor structure (without really building the factor)
			faust_unsigned_int rid; //real id
			if(this->is_transposed) {
				rid = this->size()-id-1;
			}
			else {
				rid = id;
				//dim = dim;
			}
			int nnz = L->getNbRow();
			if(id == this->size()-1)
				return nnz;
			else
			{
				nnz += this->L->getNonZeros();
				if(rid < this->size()-2)
					nnz += (this->size()-rid-1)*L->getNbRow();
			}
			return nnz;
		}

	template<typename FPP>
		bool TransformHelperPoly<FPP>::is_fact_sparse(const faust_unsigned_int id) const
		{
			return true; // all factors are sparse in a TransformHelperPoly
		}

	template<typename FPP>
		bool TransformHelperPoly<FPP>::is_fact_dense(const faust_unsigned_int id) const
		{
			return false; // all factors are sparse in a TransformHelperPoly
		}


	template<typename FPP>
		void TransformHelperPoly<FPP>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode/*=DEFAULT*/)
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}


	template<typename FPP>
		void TransformHelperPoly<FPP>::pack_factors(const faust_unsigned_int id, const PackDir dir, const int mul_order_opt_mode/*=DEFAULT*/)
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::pack_factors(const int mul_order_opt_mode/*=DEFAULT*/)
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::left(const faust_unsigned_int id, const bool copy/*=false*/) const
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::right(const faust_unsigned_int id, const bool copy/*=false*/) const
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::optimize_storage(const bool time/*=true*/)
		{
			// in general nothing better than sparse factors, just clone
			return this->clone();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::clone()
		{
			// a clone of a TransformHelperPoly is also a TransformHelperPoly (not just a TransformHelper)
			return new TransformHelperPoly<FPP>(this->size()-1, *this);
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::update_total_nnz()
		{
			//nothing to do (because of the laziness or RAII/immutability)
		}


	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::multiply(const TransformHelper<FPP, Cpu>* other) const
		{
			if(this->laziness == INSTANTIATE_COMPUTE_AND_FREE)
			{
				throw std::runtime_error("Can't multiply a FaustPoly to another Faust if highest level of lazy instantiation is enabled (INSTANTIATE_COMPUTE_AND_FREE).");
			}
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			return TransformHelper<FPP, Cpu>::multiply(other);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::multiply(FPP& scalar)
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::multiply(scalar);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelperPoly<FPP>::getNBytes() const
		{
			faust_unsigned_int nbytes = 0;
			for(int i=0;i<this->size();i++)
			{
				nbytes += MatSparse<FPP,Cpu>::getNBytes(this->get_fact_nnz(i), this->get_fact_nb_rows(i));
			}
			return nbytes;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelperPoly<FPP>::get_total_nnz() const
		{
			int total_nnz = 0;
			for(int i=0;i<this->size();i++)
				total_nnz += get_fact_nnz(i);
			return total_nnz;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::resize(faust_unsigned_int)
		{
			ERROR_ON_FAC_NUM_CHANGE();
		}

	template<typename FPP>
		std::string TransformHelperPoly<FPP>::to_string() const
		{

			stringstream str;
			str<<"Faust size ";
				str << this->get_fact_nb_rows(0) << "x" << this->get_fact_nb_cols(this->size()-1);
			auto density = (double)this->get_total_nnz()/this->getNbRow()/this->getNbCol();
			str <<", density "<< density << ", nnz_sum "<<this->get_total_nnz() << ", " << this->size() << " factor(s): "<< std::endl;
			for (int i=0 ; i<this->size() ; i++)
			{
				str << "- ";
				if(this->mul_and_combi_lin_on_gpu)
					str << "GPU ";
				str << "FACTOR " << i;
				density = (double) this->get_fact_nnz(i) / this->get_fact_nb_rows(i) / this->get_fact_nb_cols(i);
				str << Faust::MatGeneric<FPP,Cpu>::to_string(this->get_fact_nb_rows(i), this->get_fact_nb_cols(i), this->is_transposed, density, this->get_fact_nnz(i), /* is_identity */ i == this->size()-1, Sparse);
			}
			return str.str();
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelperPoly<FPP>::get_product(const int mul_order_opt_mode/*=DEFAULT*/)
		{
			//TODO: optimize: could multiply by Identity
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::get_product(mul_order_opt_mode);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::get_product(MatDense<FPP,Cpu>& prod, const int mul_order_opt_mode/*=DEFAULT*/)
		{
			//TODO: optimize: could multiply by Identity
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			TransformHelper<FPP, Cpu>::get_product(prod, mul_order_opt_mode);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
		}

	template<typename FPP>
		void TransformHelperPoly<FPP>::save_mat_file(const char* filename) const
		{
			//TODO: optimize: could instantiate only one factor at a time
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			TransformHelper<FPP, Cpu>::save_mat_file(filename);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
		}

	template<typename FPP>
		double TransformHelperPoly<FPP>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
		{
			//TODO: optimize: could instantiate only one factor at a time
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::spectralNorm(nbr_iter_max, threshold, flag);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::vertcat(const TransformHelper<FPP,Cpu>* other)
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::vertcat(other);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::horzcat(const TransformHelper<FPP,Cpu>* other)
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::horzcat(other);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		double TransformHelperPoly<FPP>::normL1() const
		{
			//TODO: optimize: could (maybe) instantiate only one factor at a time
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::normL1();
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		double TransformHelperPoly<FPP>::normFro() const
		{
			//TODO: optimize: could (maybe) instantiate only one factor at a time
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::normFro();
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		double TransformHelperPoly<FPP>::normInf() const
		{
			//TODO: optimize: could (maybe) instantiate only one factor at a time
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::normInf();
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::normalize(const int meth/*= 2 1 for 1-norm, 2 for 2-norm, MAX for inf-norm */) const
		{
			if(this->laziness == INSTANTIATE_COMPUTE_AND_FREE)
			{
				throw std::runtime_error("Can't normalize a FaustPoly if highest level of lazy instantiation is enabled (INSTANTIATE_COMPUTE_AND_FREE).");
			}
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::normalize(meth);
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::pruneout(const int nnz_tres, const int npasses/*=-1*/, const bool only_forward/*=false*/)
		{
			return this->clone(); // nothing to pruneout in a TransformHelperPoly
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::optimize_multiply(std::function<void()> f, const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/, const char* op_name/*="unamed_op"*/)
		{
			return this->clone(); // nothing to do, the multipliation is specialized in a Faust::TransformHelperPoly
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::optimize_time(const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/)
		{
			return this->clone(); // nothing to do, the multipliation is specialized in a Faust::TransformHelperPoly
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::optimize_time_full(const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/)
		{
			return this->clone(); // nothing to do, the multipliation is specialized in a Faust::TransformHelperPoly
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::optimize_time_Fv(const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/)
		{
			return this->clone(); // nothing to do, the multipliation is specialized in a Faust::TransformHelperPoly
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation/*=false*/, const bool inplace/*=false*/, const bool check_transpose/*=true*/)
		{
			if(this->laziness == INSTANTIATE_COMPUTE_AND_FREE)
			{
				throw std::runtime_error("Can't swap_cols a FaustPoly if highest level of lazy instantiation is enabled (INSTANTIATE_COMPUTE_AND_FREE).");
			}
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::swap_cols(id1, id2, permutation, inplace, check_transpose);
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperPoly<FPP>::swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation/*=false*/, const bool inplace/*=false*/, const bool check_transpose/*=true*/)
		{
			if(this->laziness == INSTANTIATE_COMPUTE_AND_FREE)
			{
				throw std::runtime_error("Can't swap_rows a FaustPoly if highest level of lazy instantiation is enabled (INSTANTIATE_COMPUTE_AND_FREE).");
			}
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelper<FPP, Cpu>::swap_rows(id1, id2, permutation, inplace, check_transpose);
			return ret;
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperPoly<FPP>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			auto this_ = const_cast<TransformHelperPoly<FPP>*>(this);
			this_->basisChebyshev_fact_all();
			auto ret = TransformHelperGen<FPP, Cpu>::slice(start_row_id, end_row_id, start_col_id, end_col_id);
			if(this_->laziness == INSTANTIATE_COMPUTE_AND_FREE)
				this_->basisChebyshev_free_fact_all();
			return ret;
		}

}
