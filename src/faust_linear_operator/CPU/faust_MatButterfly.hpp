
namespace Faust
{
	template<typename FPP>
	MatButterfly<FPP, Cpu>::MatButterfly(const MatSparse<FPP, Cpu> &factor, int level)
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
	MatButterfly<FPP, Cpu>::MatButterfly(const MatButterfly& src)
	{
		*this = src;
	}


	template<typename FPP>
	MatButterfly<FPP, Cpu>& MatButterfly<FPP, Cpu>::operator=(const MatButterfly& src)
	{
		D1 = src.D1;
		D2 = src.D2;
		D2T = src.D2T;
		subdiag_ids = src.subdiag_ids;
#ifdef USE_PYTHONIC
		subdiag_ids_ptr = new long[src.getNbRow()];
		copy(src.subdiag_ids_ptr, src.subdiag_ids_ptr + src.getNbRow(), subdiag_ids_ptr);
#endif
		level = src.level;
		return *this;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::Display() const
	{
		//TODO: adjust consistently with MatGeneric Display (using to_string)
		std::cout << "MatButterfly on CPU: ";
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
		void MatButterfly<FPP, Cpu>::init_transpose()
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
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::Clone(const bool isOptimize) const
	{
		return new MatButterfly<FPP, Cpu>(*this);
	}

	template<typename FPP>
	Vect<FPP, Cpu> MatButterfly<FPP, Cpu>::multiply(const Vect<FPP,Cpu> & x) const
	{
		Vect<FPP, Cpu> z(x.size());
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), false); //TODO: conjugate
		return z;
	}


	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(Vect<FPP,Cpu> & x, char opThis) const
	{
		Vect<FPP, Cpu> z(x.size());
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), opThis != 'N'); //TODO: conjugate
		x = z;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
	{

		MatDense<FPP, Cpu> Y(this->getNbRow(), M.getNbCol());
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(M.getData(), M.getNbCol(), Y.getData(), Y.getNbRow(), 'T' == opThis);
		M = Y;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(MatSparse<FPP, Cpu>& M, char opThis) const
	{
		MatDense<FPP, Cpu> Y(M);
		this->multiply(Y, opThis);
		M = Y;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiplyRight(MatSparse<FPP, Cpu> const& M)
	{
		throw std::runtime_error("multiplyRight is not supported on a MatButterfly because the matrix must stay a butterfly matrix.");
	}


	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(const FPP* x, FPP* y, size_t size, bool transpose/* = false*/)
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
		void  MatButterfly<FPP, Cpu>::multiply(const FPP* X, int X_ncols, FPP* Y, size_t Y_nrows, bool transpose/* = false*/)
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
		void MatButterfly<FPP, Cpu>::faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const
		{
			//TODO
		}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::transpose()
	{
		//TODO
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::conjugate(const bool eval)
	{
		//TODO
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::adjoint()
	{
		//TODO
	}

	template<typename FPP>
	faust_unsigned_int MatButterfly<FPP, Cpu>::getNonZeros()const
	{
		//TODO
	}

	template<typename FPP>
	size_t MatButterfly<FPP, Cpu>::getNBytes() const
	{
		//TODO
	}

	template<typename FPP>
	MatType MatButterfly<FPP, Cpu>::getType() const
	{
		//TODO
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::operator*=(const FPP alpha)
	{
		//TODO
	}

	template<typename FPP>
	matvar_t* MatButterfly<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate, const char *var_name) const
	{
		//TODO
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::normL1(const bool transpose) const
	{
		//TODO
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::norm() const
	{
		//TODO
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose) const
	{
		//TODO
	}

	template<typename FPP>
	Vect<FPP,Cpu> MatButterfly<FPP, Cpu>::get_col(faust_unsigned_int id) const
	{
		//TODO
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
	{
		//TODO
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
	{
		//TODO
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
	{
		//TODO
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
	{
		//TODO
	}

	template<typename FPP>
	std::list<std::pair<int,int>> MatButterfly<FPP, Cpu>::nonzeros_indices() const
	{
		//TODO
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::setZeros()
	{
		//TODO
	}

	template<typename FPP>
	bool MatButterfly<FPP, Cpu>::containsNaN() const
	{
		//TODO
	}

	template<typename FPP>
	const FPP& MatButterfly<FPP, Cpu>::operator()(faust_unsigned_int i, faust_unsigned_int j) const
	{
		//TODO
	}

}
