#include <cmath>
namespace Faust
{

	template<typename FPP>
		TransformHelperButterfly<FPP,Cpu>::TransformHelperButterfly(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ /*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/) : TransformHelper<FPP, Cpu>(facts, lambda_, optimizedCopy, cloning_fact, internal_call)
	{
		int i = 0;
		auto size = this->getNbRow();
//		for(auto csr_fac: facts)
		// use rather recorded factors in the Faust::Transform because one might have been multiplied with lambda_
		auto log2nf = 1 << (this->size() - 1);
		has_permutation = (log2nf - this->getNbRow()) == 0;
		auto end_it = has_permutation?this->end()-1:this->end();
		for(auto csr_fac_it = this->begin(); csr_fac_it != end_it; csr_fac_it++)
		{
			auto csr_fac = *csr_fac_it;
			opt_factors.insert(opt_factors.begin(), std::make_shared<MatButterfly<FPP, Cpu>>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac), i++));
		}
		if(has_permutation)
		{
			// set the permutation factor
			auto csr_fac = *(this->end()-1);
			perm = MatPerm<FPP, Cpu>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac));
		}
	}


	template<typename FPP>
		TransformHelperButterfly<FPP,Cpu>::TransformHelperButterfly(const TransformHelperButterfly<FPP,Cpu>* th, bool transpose, bool conjugate) : TransformHelper<FPP,Cpu>()
	{
		// don't use parent transpose/conjugate ctor because it fails with other variadic template ctor taking over
		this->transform = th->transform;
		this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
		this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		this->copy_slice_state(*th, transpose);
		this->copy_fancy_idx_state(*th, transpose);
		this->copy_mul_mode_state(*th);
		opt_factors = th->opt_factors;
		if(this->is_transposed)
			// factors are reversed if the faust is transposed
			std::reverse(opt_factors.begin(), opt_factors.end());

		if(th->has_permutation)
		{
			// set the transpose permutation factor
			// TODO: refactor with previous ctor
			auto csr_fac = *(this->end()-1);
			auto size = csr_fac->getNbRow();
			MatSparse<FPP, Cpu> sp_fac = MatSparse<FPP, Cpu>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac));
			perm = MatPerm<FPP, Cpu>(sp_fac);
			if(transpose)
				perm.transpose();
			if(conjugate)
				perm.conjugate();
		}
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelperButterfly<FPP,Cpu>::transpose()
	{
		return new TransformHelperButterfly<FPP,Cpu>(this, true, false);
	}

	template<typename FPP>
		std::string TransformHelperButterfly<FPP,Cpu>::to_string() const
		{
			auto str = TransformHelper<FPP,Cpu>::to_string();
			str += "- Butterfly structure optimized\r\n";
			return str;
		}


	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelperButterfly<FPP,Cpu>::fourierFaust(unsigned int n, const bool norma)
		{

			std::vector<MatGeneric<FPP,Cpu>*> factors(n+1);
			TransformHelper<FPP,Cpu>* fourierFaust = nullptr;
			try
			{
				fft_factors(n, factors);
				FPP alpha = norma?FPP(1/sqrt((double)(1 << n))):FPP(1.0);
				fourierFaust = new TransformHelperButterfly<FPP, Cpu>(factors, alpha, false, false, /* internal call */ true);
			}
			catch(std::bad_alloc e)
			{
				//nothing to do, out of memory, return nullptr
			}
			return fourierFaust;
		}


	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelperButterfly<FPP,Cpu>::optFaust(const TransformHelper<FPP, Cpu>* F)
		{
			//TODO: verify a few assertions to detect if F factors do not match a butterfly structure
			std::vector<MatGeneric<FPP,Cpu>*> factors(F->size());
			copy(F->begin(), F->end(), factors.begin());
			auto oF = new TransformHelperButterfly<FPP, Cpu>(factors, FPP(1.0), false, false, /* internal call */ true);
			return oF;
		}


	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperButterfly<FPP,Cpu>::multiply(const Vect<FPP, Cpu>& x)
		{
			Vect<FPP, Cpu> y(this->getNbRow());
			multiply(x.getData(), y.getData());
			return y;
		}

	template<typename FPP>
		void TransformHelperButterfly<FPP,Cpu>::multiply(const FPP* x, FPP* y)
		{
//			std::cout << "THButterfly multiply vec" << std::endl;
//			std::cout << "transpose:" << this->is_transposed << " has_permutation:" << has_permutation << std::endl;
			auto size = this->getNbRow();
			VecMap x_vec(const_cast<FPP*>(x), size); // const_cast is harmless
			Vect<FPP, Cpu> z(size);
			int i = 0;
			if(has_permutation && ! this->is_transposed)
			{
				perm.multiply(x, y, size, false, this->is_conjugate);
			}
			else
			{
				MatButterfly<FPP, Cpu>& fac = * opt_factors[0];
				fac.multiply(x, z.getData(), this->getNbRow(), this->is_transposed);
				i = 1;
			}

			while(i < opt_factors.size())
			{
				MatButterfly<FPP, Cpu>& fac = * opt_factors[i];
				if(i & 1)
					fac.multiply(z.getData(), y, this->getNbRow(), this->is_transposed);
				else
					fac.multiply(y, z.getData(), this->getNbRow(), this->is_transposed);
				i++;
			}
			if(has_permutation && this->is_transposed)
			{
				if(i & 1)
				{
					// the vector is in z
					perm.multiply(z.getData(), y, this->getNbRow(), true, this->is_conjugate);
				}
				else
				{
					// the vector is in y
					perm.multiply(y, z.getData(), this->getNbRow(), true, this->is_conjugate);
					memcpy(y, z.getData(), size*sizeof(FPP));
				}
			}
			else
				if(i & 1)
					memcpy(y, z.getData(), size*sizeof(FPP));
		}

	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperButterfly<FPP,Cpu>::multiply(const FPP* x)
		{
			Vect<FPP, Cpu> y(this->getNbRow());
			multiply(x, y.getData());
			return y;
		}

	template<typename FPP>
			void  TransformHelperButterfly<FPP,Cpu>::multiply(const FPP* X, int X_ncols, FPP* Y)
			{
				using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
				MatMap X_mat(const_cast<FPP*>(X) /* harmless, no modification*/, this->getNbCol(), X_ncols);
				auto Z = new FPP[this->getNbRow()*X_ncols];
				int i = 0;
				MatMap Y_mat(Y, this->getNbRow(), X_ncols);
				if(has_permutation && ! this->is_transposed)
				{
					perm.multiply(X, X_ncols, Y, this->getNbRow(), false, this->is_conjugate);
				}
				else
				{
					MatButterfly<FPP, Cpu>& fac = * opt_factors[0];
					fac.multiply(X, X_mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i = 1;
				}

				while(i < opt_factors.size())
				{
					MatButterfly<FPP, Cpu>& fac = const_cast<MatButterfly<FPP, Cpu>&>(* opt_factors[i]);

					if(i & 1)
						fac.multiply(Z, Y_mat.cols(), Y, this->getNbRow(), this->is_transposed);
					else
						fac.multiply(Y, Y_mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i++;
				}

				if(has_permutation && this->is_transposed)
				{
					if(i & 1)
						perm.multiply(Z, X_ncols, Y, this->getNbRow(), true, this->is_conjugate);
					else
					{
						perm.multiply(Y, X_ncols, Z, this->getNbRow(), true, this->is_conjugate);
						memcpy(Y, Z, sizeof(FPP)*this->getNbRow()*X_ncols);
					}
				}
				else
					if(i & 1)
						memcpy(Y, Z, sizeof(FPP)*this->getNbRow()*X_ncols);
				delete[] Z;
			}

	template<typename FPP>
			MatDense<FPP, Cpu>  TransformHelperButterfly<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> &X)
			{
				MatDense<FPP, Cpu> Y(this->getNbRow(), X.getNbCol());
				multiply(X.getData(), X.getNbCol(), Y.getData());
				return Y;
			}

	template<typename FPP>
			MatDense<FPP, Cpu>  TransformHelperButterfly<FPP,Cpu>::multiply(const MatSparse<FPP,Cpu> &X)
			{
				MatDense<FPP, Cpu> Y(this->getNbRow(), X.getNbCol());
				auto Z = new FPP[this->getNbRow()*X.getNbCol()];
				int i = 0;
				if(has_permutation && ! this->is_transposed)
				{
					//TODO: conjugate
					Y = X;
					perm.multiply(Y, 'N');
				}
				else
				{
					Y = X;
				}
#ifdef FAUST_BUTTERFLY_MULTILPY_SINGLE_BUFFER
				// consume less memory but with more copies
				//					for(auto fac: opt_factors)
				while(i < opt_factors.size())
				{
					MatButterfly<FPP, Cpu>& fac = * opt_factors[i];
					Y = fac.multiply(Y, this->is_transposed);
					i++;
				}
#else
				// double buffering consumes more memory but spare copies (faster if enough memory)
				//TODO: factorize with MatDense product
				while(i < opt_factors.size())
				{
//					ButterflyMat<FPP, Cpu>& fac = * opt_factors[i];
					MatButterfly<FPP, Cpu>& fac = * opt_factors[i];
					if(i & 1)
						fac.multiply(Z, Y.mat.cols(), Y.getData(), this->getNbRow(), this->is_transposed);
					else
						fac.multiply(Y.getData(), Y.mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i++;
				}
				if(has_permutation && this->is_transposed)
				{
					if(i & 1)
						perm.multiply(Z, X.getNbCol(), Y.getData(), this->getNbRow(), true, this->is_conjugate);
					else
					{
						perm.multiply(Y.getData(), X.getNbCol(), Z, this->getNbRow(), true, this->is_conjugate);
						memcpy(Y.getData(), Z, sizeof(FPP)*this->getNbRow()*X.getNbCol());
					}
				}
				else if(i & 1)
					memcpy(Y.getData(), Z, sizeof(FPP)*this->getNbRow()*X.getNbCol());
				delete[] Z;
#endif
				return Y;
			}

}
