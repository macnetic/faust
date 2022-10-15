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
			opt_factors.insert(opt_factors.begin(), std::make_shared<ButterflyMat<FPP, Cpu>>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac), i++));
		}
		if(has_permutation)
		{
			// set the permutation factor
			auto csr_fac = *(this->end()-1);
			D.resize(size);
			perm_d_ptr = D.diagonal().data();
			// only a setOnes should be enough because this is a permutation matrix (but it could be normalized)
			memcpy(perm_d_ptr, dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac)->getValuePtr(), size*sizeof(FPP));
//			auto bitrev_perm_ids = new unsigned int[size];
//			iota(bitrev_perm_ids, bitrev_perm_ids+size, 0);
//			bit_rev_permu(facts.size()-1, bitrev_perm_ids);
			bitrev_perm.resize(size);
//			copy(bitrev_perm_ids, bitrev_perm_ids+size, bitrev_perm.begin());
//			delete[] bitrev_perm_ids;
			//TODO: rename bitrev_perm to something generic (this class is not limited to the FFT)
			copy(dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac)->getColInd(), dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac)->getColInd()+size, bitrev_perm.begin());
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
			if(transpose)
				sp_fac.transpose();
			if(conjugate)
				sp_fac.conjugate();
			D.resize(size);
			perm_d_ptr = D.diagonal().data();
			memcpy(perm_d_ptr, sp_fac.getValuePtr(), size*sizeof(FPP));
			bitrev_perm.resize(size);
			copy(sp_fac.getColInd(), sp_fac.getColInd()+size, bitrev_perm.begin());
			this->has_permutation = true;
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
				VecMap y_vec(y, size);
				if(x == y)
				{
					// an intermediate vector is needed to index x
					auto x_ = new FPP[size];
#pragma omp parallel for
					for(int i=0;i < size; i++)
						x_[i] = x[bitrev_perm[i]];
#pragma omp parallel for
					for(int i=0;i < size; i++)
						y[i] = perm_d_ptr[i] * x_[i];
					delete[] x_;
				}
				else
#define BUTTERFLY_MUL_VEC_OMP_LOOP
#ifdef BUTTERFLY_MUL_VEC_OMP_LOOP
					// faster
#pragma omp parallel for
					for(int i=0;i < this->getNbRow(); i++)
						y[i] = perm_d_ptr[i] * x[bitrev_perm[i]];
#else
				y_vec = D * x_vec(bitrev_perm);
#endif
			}
			else
			{
				ButterflyMat<FPP, Cpu>& fac = * opt_factors[0];
				fac.multiply(x, z.getData(), this->getNbRow(), this->is_transposed);
				i = 1;
			}

			while(i < opt_factors.size())
			{
				ButterflyMat<FPP, Cpu>& fac = * opt_factors[i];
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
					//TODO: refactor eltwise mul
                   #pragma omp parallel for
					for(int i=0;i < this->getNbRow(); i++)
						y[i] = perm_d_ptr[i] * z[bitrev_perm[i]];
				}
				else
				{
					// the vector is in y
                   #pragma omp parallel for
					for(int i=0;i < this->getNbRow(); i++)
						z[i] = perm_d_ptr[i] * y[bitrev_perm[i]];
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
				if(has_permutation)
				{
#if defined(BUTTERFLY_MUL_MAT_OMP_LOOP) || ! EIGEN_VERSION_AT_LEAST(3, 4, 0)
					// this is slower
#pragma parallel omp for
					for(int i=0;i < this->getNbRow(); i ++)
						Y_mat.row(i) = X_mat.row(bitrev_perm[i]) * perm_d_ptr[i];
#else
					Y_mat = D * X_mat(bitrev_perm, Eigen::placeholders::all);
#endif
				}
				else
				{
					ButterflyMat<FPP, Cpu>& fac = * opt_factors[0];
					fac.multiply(X, X_mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i = 1;
				}

				while(i < opt_factors.size())
					//			for(auto fac: opt_factors)
				{
					ButterflyMat<FPP, Cpu>& fac = * opt_factors[i];

					if(i & 1)
						fac.multiply(Z, Y_mat.cols(), Y, this->getNbRow(), this->is_transposed);
					else
						fac.multiply(Y, Y_mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i++;
				}
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
				if(has_permutation)
				{
#pragma omp parallel for
					for(int i=0;i < this->getNbRow(); i ++)
						Y.mat.row(i) = X.mat.row(bitrev_perm[i]) * perm_d_ptr[i];

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
					ButterflyMat<FPP, Cpu>& fac = * opt_factors[i];
					Y = fac.multiply(Y, this->is_transposed);
					i++;
				}
#else
				// double buffering consumes more memory but spare copies (faster if enough memory)
				//TODO: factorize with MatDense product
				while(i < opt_factors.size())
				{
					ButterflyMat<FPP, Cpu>& fac = * opt_factors[i];
					if(i & 1)
						fac.multiply(Z, Y.mat.cols(), Y.getData(), this->getNbRow(), this->is_transposed);
					else
						fac.multiply(Y.getData(), Y.mat.cols(), Z, this->getNbRow(), this->is_transposed);
					i++;
				}

				if(i & 1)
					memcpy(Y.getData(), Z, sizeof(FPP)*this->getNbRow()*X.getNbCol());
				delete[] Z;
#endif
				return Y;
			}

}


namespace Faust
{
	template<typename FPP>
	ButterflyMat<FPP, Cpu>::ButterflyMat(const MatSparse<FPP, Cpu> &factor, int level)
	{
		// build a d1, d2 pair from the butterfly factor
		auto size = factor.getNbRow();
		D1.resize(size);
		D2.resize(size);
		D2T.resize(0); // not used defaultly
		FPP *d1_ptr, *d2_ptr;
		d1_ptr = D1.diagonal().data();
		d2_ptr = D2.diagonal().data();
		auto d_offset = size >> (level+1);
		auto data = factor.getValuePtr();
		auto rowptr = factor.getRowPtr();
		for(int i=0;i < size; i++)
		{
			if((i / d_offset) & 1)
			{
				// d2 coeff is the first elt of row i
				d2_ptr[i] = data[rowptr[i]];
				d1_ptr[i] = data[rowptr[i]+1]; // diag elt is just after
			}
			else
			{
				// d2 coeff is the last elt of row i
				d2_ptr[i] = data[rowptr[i+1]-1];
				d1_ptr[i] = data[rowptr[i+1]-2]; // diag elt is just before
			}
		}
		std::vector<int> seq(size);
		iota(seq.begin(), seq.end(), 0);
		subdiag_ids.resize(size);
		for(int i = 0;i < size; i += d_offset * 2)
		{
			copy(seq.begin()+i+d_offset, seq.begin()+i+2*d_offset, subdiag_ids.begin()+i);
			copy(seq.begin()+i, seq.begin()+i+d_offset, subdiag_ids.begin()+i+d_offset);
		}
#ifdef USE_PYTHONIC
		subdiag_ids_ptr = new long[size];
		copy(subdiag_ids.begin(), subdiag_ids.end(),subdiag_ids_ptr);
#endif
		this->level = level;
	}

	template<typename FPP>
		void ButterflyMat<FPP, Cpu>::init_transpose()
		{
			//TODO: simplify in case of symmetric matrix (it happens for the FFT)
			auto size = D2.rows();
			if(D2T.size() != 0)
				// already done
				return;
			FPP *d1_ptr, *d2_ptr, *d2t_ptr;
			d1_ptr = D1.diagonal().data();
			d2_ptr = D2.diagonal().data();
			D2T.resize(size);
			d2t_ptr = D2T.diagonal().data();

			auto d_offset = size >> (level+1);
			// D1 doesn't change
			// swap every pair of D2 contiguous blocks to form D2T
			for(int i = 0;i < size; i += d_offset * 2)
			{
				// swap two next blocks of size d_offset into d2t_ptr
				copy(d2_ptr + i, d2_ptr + i + d_offset, d2t_ptr + i + d_offset);
				copy(d2_ptr + i + d_offset, d2_ptr + i + 2 * d_offset, d2t_ptr + i);
			}
		}

	template<typename FPP>
	void ButterflyMat<FPP, Cpu>::Display() const
	{
		std::cout << "ButterflyMat on CPU: ";
		std::cout << "D1: ";
		std::cout << D1.diagonal() << std::endl;
		std::cout << "D2: ";
		std::cout << D2.diagonal() << std::endl;
		if(D2T.size() > 0)
		{
			std::cout << "D2T: ";
			std::cout << D2T.diagonal() << std::endl;
		}
		cout << "subdiag_ids: ";
		for(int i=0;i < subdiag_ids.size();i++)
			cout << subdiag_ids[i] << " ";
		cout << std::endl;
	}

	template<typename FPP>
	Vect<FPP, Cpu> ButterflyMat<FPP, Cpu>::multiply(const Vect<FPP, Cpu>& x, bool transpose/* = false*/)
	{
		Vect<FPP, Cpu> z(x.size());
		multiply(x.getData(), z.getData(), x.size(), transpose);
		return z;
	}

	template<typename FPP>
	void ButterflyMat<FPP, Cpu>::multiply(const FPP* x, FPP* y, size_t size, bool transpose/* = false*/)
	{
		DiagMat &D2 = this->D2;
		const FPP *d1_ptr, *d2_ptr;
		d1_ptr = D1.diagonal().data();
		if(transpose)
		{
			init_transpose(); // no cost if already initialized
			d2_ptr = D2T.diagonal().data();
			D2 = D2T;
		}
		else
			d2_ptr = D2.diagonal().data();
#ifdef USE_PYTHONIC
		auto xArray = arrayFromBuf1D(x, size);
		auto d1Array = arrayFromBuf1D(d1_ptr, size);
		auto d2Array = arrayFromBuf1D(d2_ptr, size);
		auto x_ids = arrayFromBuf1D(subdiag_ids_ptr, size);
//		auto yArray = __pythran_ButFactor_matmul::__matmul__()(d1Array, d2Array, xArray, x_ids);
		auto yArray = arrayFromBuf1D(y, size);
		yArray = pythonic::operator_::add(pythonic::operator_::mul(d1Array, xArray), pythonic::operator_::mul(d2Array, xArray[x_ids]));
		memcpy(y, yArray.buffer, sizeof(FPP)*size);
#endif
#define BMAT_MULTIPLY_VEC_OMP_LOOP
#ifdef BMAT_MULTIPLY_VEC_OMP_LOOP
		#pragma omp parallel for
		for(int i=0;i < size; i++)
			y[i] = d1_ptr[i] * x[i] + d2_ptr[i] * x[subdiag_ids[i]];
#else
		// this is slower
		VecMap x_vec(const_cast<FPP*>(x), size); // const_cast is harmless
		VecMap y_vec(y, size); // const_cast is harmless
		y_vec = D1 * x_vec + D2 * x_vec(subdiag_ids);
#endif
	}

	template<typename FPP>
		void  ButterflyMat<FPP, Cpu>::multiply(const FPP* X, int X_ncols, FPP* Y, size_t Y_nrows, bool transpose/* = false*/)
		{
			using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
			MatMap X_mat(const_cast<FPP*>(X) /* harmless, no modification*/, Y_nrows, X_ncols);
			MatMap Y_mat(Y, Y_nrows, X_ncols);
			// --------- TODO refactor this block with overload of multiply
			DiagMat &D2 = this->D2;
			const FPP *d1_ptr, *d2_ptr;
			d1_ptr = D1.diagonal().data();
			if(transpose)
			{
				init_transpose(); // no cost if already initialized
				d2_ptr = D2T.diagonal().data();
				D2 = D2T;
			}
			else
				d2_ptr = D2.diagonal().data();
			// ---------------
#if defined(BUTTERFLY_MUL_MAT_OMP_LOOP) || ! EIGEN_VERSION_AT_LEAST(3, 4, 0)
			// this is slower
			#pragma omp parallel for
			for(int i=0;i < Y_nrows; i++)
				Y_mat.row(i) = d1_ptr[i] * X_mat.row(i) + d2_ptr[i] * X_mat.row(subdiag_ids[i]);
#else
			Y_mat = D1 * X_mat + D2 * X_mat(subdiag_ids, Eigen::placeholders::all);
#endif
		}

	template<typename FPP>
		MatDense<FPP, Cpu> ButterflyMat<FPP, Cpu>::multiply(const MatDense<FPP,Cpu> &X, bool transpose/* = false*/)
		{
			MatDense<FPP, Cpu> Y(X.getNbRow(), X.getNbCol());
			multiply(X.getData(), X.getNbCol(), Y.getData(), X.getNbRow(), transpose);
			return Y;
		}
}
