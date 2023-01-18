#include "faust_conj.h"
namespace Faust
{

	//TODO: dtor?
	//
	template<typename FPP>
	MatPerm<FPP, Cpu>::MatPerm()
	{
#ifdef USE_PYTHONIC
		perm_ids_ptr = nullptr;
#endif
		is_transp = false;
	}

	template<typename FPP>
	MatPerm<FPP, Cpu>::MatPerm(const MatSparse<FPP, Cpu> &sp_mat) : MatPerm()
	{
		//TODO: test sp_mat is a permutation matrix (at least nrows == ncols)
		auto size = sp_mat.getNbRow();
		D.resize(size);

		FPP *perm_d_ptr = D.diagonal().data();

		// only a setOnes should be enough because this is a permutation matrix (but it could be normalized)
		memcpy(perm_d_ptr, sp_mat.getValuePtr(), size*sizeof(FPP));

		perm_ids.resize(size);
		std::copy(sp_mat.getColInd(), sp_mat.getColInd()+size, perm_ids.begin());
	}

	template<typename FPP>
	MatPerm<FPP, Cpu>::MatPerm(const MatPerm& src) : MatPerm()
	{
		*this = src;
	}


	template<typename FPP>
	MatPerm<FPP, Cpu>& MatPerm<FPP, Cpu>::operator=(const MatPerm& src)
	{
		D = src.D;
		perm_ids = src.perm_ids;
#ifdef USE_PYTHONIC
		perm_ids_ptr = new long[src.getNbRow()];
		std::copy(src.perm_ids_ptr, src.perm_ids_ptr + src.getNbRow(), perm_ids_ptr);
#endif
		if(src.DT.size() > 0)
		{
			DT = src.DT;
			perm_ids_T.resize(src.perm_ids_T.size());
			std::copy(src.perm_ids_T.begin(), src.perm_ids_T.end(), perm_ids_T.begin());
		}
		is_transp = src.is_transp;
		return *this;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::Display() const
	{
		//TODO: adjust consistently with MatGeneric Display (using to_string)
		std::cout << "MatPerm on CPU: ";
		std::cout << "D: ";
		std::cout << D.diagonal() << std::endl;
		std::cout << "perm_ids: ";
		for(int i=0;i < perm_ids.size();i++)
			std::cout << perm_ids[i] << " ";
		std::cout << std::endl;
		if(DT.size() > 0)
		{
			std::cout << "DT: ";
			std::cout << DT.diagonal() << std::endl;
			std::cout << "perm_ids_T: ";
			for(int i=0;i < perm_ids.size();i++)
				std::cout << perm_ids_T[i] << " ";
			std::cout << std::endl;
		}
	}

	template<typename FPP>
		void MatPerm<FPP, Cpu>::init_transpose()
		{
			//TODO: simplify in case of symmetric matrix (it happens for the FFT)
			//TODO: optimize if it's a proper permutation (only ones)
			if(DT.size() == 0)
			{
				DT.resize(D.rows());
				perm_ids_T.resize(perm_ids.size());
				auto dt_ptr = DT.diagonal().data();
				auto d_ptr = D.diagonal().data();
				for(int i=0;i < D.rows(); i++)
				{
					auto cid = perm_ids[i];
					perm_ids_T[cid] = i;
					dt_ptr[cid] = d_ptr[i];
				}
			}
		}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatPerm<FPP, Cpu>::Clone(const bool isOptimize) const
	{
		return new MatPerm<FPP, Cpu>(*this);
	}

	template<typename FPP>
	Vect<FPP, Cpu> MatPerm<FPP, Cpu>::multiply(const Vect<FPP,Cpu> & x) const
	{
		Vect<FPP, Cpu> z(x.size());
		const_cast<MatPerm<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), false); //TODO: conjugate
		return z;
	}


	template<typename FPP>
	void MatPerm<FPP, Cpu>::multiply(Vect<FPP,Cpu> & x, char opThis) const
	{
		Vect<FPP, Cpu> z(x.size());
		const_cast<MatPerm<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), opThis != 'N'); //TODO: conjugate
		x = z;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
	{
		MatDense<FPP, Cpu> Y(this->getNbRow(), M.getNbCol());
		const_cast<MatPerm<FPP, Cpu>*>(this)->multiply(M.getData(), M.getNbCol(), Y.getData(), Y.getNbRow(), 'N' != opThis, 'H' == opThis || 'C' == opThis);
		M = Y;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::multiply(MatSparse<FPP, Cpu>& M, char opThis) const
	{
		MatDense<FPP, Cpu> Y(M);
		this->multiply(Y, opThis);
		M = Y;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::multiplyRight(MatSparse<FPP, Cpu> const& M)
	{
		throw std::runtime_error("multiplyRight is not supported on a MatPerm because the matrix must stay a butterfly matrix.");
	}


	template<typename FPP>
	void MatPerm<FPP, Cpu>::multiply(const FPP* x, FPP* y, size_t size, bool transpose/* = false*/, bool conjugate/* = false*/)
	{
		DiagMat &D = this->D;
		const FPP *d_ptr;
		auto do_transp = transpose ^ is_transp;
		if(do_transp)
		{
			init_transpose(); // no cost if already initialized
			d_ptr = DT.diagonal().data();
		}
		else
			d_ptr = D.diagonal().data();
#ifdef USE_PYTHONIC
		auto xArray = arrayFromBuf1D(x, size);
		auto d2Array = arrayFromBuf1D(d_ptr, size);
		auto x_ids = arrayFromBuf1D(perm_ids_ptr, size);
		//		auto yArray = __pythran_ButFactor_matmul::__matmul__()(d1Array, d2Array, xArray, x_ids);
		auto yArray = arrayFromBuf1D(y, size);
		yArray = pythonic::operator_::mul(d2Array, xArray[x_ids]);
		memcpy(y, yArray.buffer, sizeof(FPP)*size);
#endif
#define BMAT_MULTIPLY_VEC_OMP_LOOP
#ifdef BMAT_MULTIPLY_VEC_OMP_LOOP
#pragma omp parallel for
		for(int i=0;i < size; i++)
			if(conjugate)
				y[i] = Faust::conj(d_ptr[i]) * x[perm_ids[i]];
			else
				y[i] = d_ptr[i] * x[perm_ids[i]];
#else
		// this is slower
		VecMap x_vec(const_cast<FPP*>(x), size); // const_cast is harmless
		VecMap y_vec(y, size); // const_cast is harmless
		if(do_transp)
			if(conjugate)
				y_vec = diag_conj(DT) * x_vec(perm_ids);
			else
				y_vec = DT * x_vec(perm_ids);
		else
			if(conjugate)
				y_vec = diag_conj(D) * x_vec(perm_ids);
			else
				y_vec = D * x_vec(perm_ids);
#endif
	}

	template<typename FPP>
		void  MatPerm<FPP, Cpu>::multiply(const FPP* X, int X_ncols, FPP* Y, size_t Y_nrows, bool transpose/* = false*/, bool conjugate /* = false*/)
		{
			using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
			MatMap X_mat(const_cast<FPP*>(X) /* harmless, no modification*/, Y_nrows, X_ncols);
			MatMap Y_mat(Y, Y_nrows, X_ncols);
			// --------- TODO refactor this block with overload of multiply
			const FPP *d_ptr;
			auto do_transp = transpose ^ is_transp;
			if(do_transp)
			{
				init_transpose(); // no cost if already initialized
				d_ptr = DT.diagonal().data();
			}
			else
				d_ptr = D.diagonal().data();
			// ---------------
#if defined(BUTTERFLY_MUL_MAT_OMP_LOOP) || ! EIGEN_VERSION_AT_LEAST(3, 4, 0)
			// this is slower
			#pragma omp parallel for
			for(int i=0;i < Y_nrows; i++)
				if(conjugate)
					Y_mat.row(i) = Faust::conj(d_ptr[i]) * X_mat.row(perm_ids[i]);
				else
					Y_mat.row(i) = d_ptr[i] * X_mat.row(perm_ids[i]);
#else
			// TODO: refactor the exp with a macro diag_prod
			if(do_transp)
				if(conjugate)
					Y_mat = diag_conj(DT) * X_mat(perm_ids, Eigen::placeholders::all);
				else
					Y_mat = DT * X_mat(perm_ids, Eigen::placeholders::all);
			else
				if(conjugate)
					Y_mat = diag_conj(D) * X_mat(perm_ids, Eigen::placeholders::all);
				else
					Y_mat = D * X_mat(perm_ids, Eigen::placeholders::all);
#endif
		}


	template<typename FPP>
		void MatPerm<FPP, Cpu>::faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const
		{
			std::runtime_error op_except("Invalid operation letter in MatPerm::faust_gemm; it must be 'N', 'T' or 'H'.");
			if(beta == FPP(0))
			{
				if(typeB == 'N')
				{
					C = B;
					multiply(C, typeA);
				}
				else if(typeB == 'T')
				{
					auto C = B;
					C.transpose();
					multiply(C, typeA);
				}
				else if(typeB == 'H')
				{
					auto C = B;
					C.adjoint();
					multiply(C, typeA);
				}
				else
					throw op_except;
				C *= alpha;
			}
			else // beta != 0
			{
				C *= beta;
				MatDense<FPP, Cpu> Bc(B); //copy
				if(alpha != FPP(0))
					Bc *= alpha;
				if(typeB == 'N')
				{
					multiply(Bc, typeA);
				}
				else if(typeB == 'T')
				{
					Bc.transpose();
					multiply(Bc, typeA);
				}
				else if(typeB == 'H')
				{

					Bc.adjoint();
					multiply(Bc, typeA);
				}
				else
					throw op_except;
				C += Bc;
			}
		}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::transpose()
	{
		init_transpose(); // free cost if already called once
		is_transp = ! is_transp;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::conjugate(const bool eval)
	{
		auto size = getNbRow();
		VecMap d_vec(const_cast<FPP*>(D.diagonal().data()), size); // const_cast is harmless
		d_vec = d_vec.conjugate();
		if(DT.size() > 0)
		{
			VecMap d_vec(const_cast<FPP*>(DT.diagonal().data()), size);
			d_vec = d_vec.conjugate();
		}
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::adjoint()
	{
		transpose();
		conjugate();
	}

	template<typename FPP>
	faust_unsigned_int MatPerm<FPP, Cpu>::getNonZeros()const
	{
		return getNbRow();
	}

	template<typename FPP>
	size_t MatPerm<FPP, Cpu>::getNBytes() const
	{
		return (D.rows() + (DT.size() != 0?DT.rows():0)) * sizeof(FPP) + perm_ids.size() * sizeof(int);
	}

	template<typename FPP>
	MatType MatPerm<FPP, Cpu>::getType() const
	{
		return Perm;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::operator*=(const FPP alpha)
	{
		D = D * alpha;
		if(DT.size() != 0)
			DT = DT * alpha;
	}

	template<typename FPP>
	matvar_t* MatPerm<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate, const char *var_name) const
	{
		throw std::runtime_error("Saving MatPerm to a .mat file is not supported.");
	}

	template<typename FPP>
	Real<FPP> MatPerm<FPP, Cpu>::normL1(const bool transpose/*=false*/) const
	{
		return toMatSparse().normL1(transpose);
	}

	template<typename FPP>
	Real<FPP> MatPerm<FPP, Cpu>::norm() const
	{
		auto s = (D.diagonal().array() * D.diagonal().conjugate().array()).sum();
		return Faust::fabs(std::sqrt(s));
	}

	template<typename FPP>
	Real<FPP> MatPerm<FPP, Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose) const
	{
		return toMatSparse().normL1(col_id, transpose);
	}

	template<typename FPP>
	Vect<FPP,Cpu> MatPerm<FPP, Cpu>::get_col(faust_unsigned_int id) const
	{
		return toMatSparse().get_col(id);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatPerm<FPP, Cpu>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
	{
		return toMatSparse().get_cols(col_id_start, num_cols);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatPerm<FPP, Cpu>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
	{
		return toMatSparse().get_rows(row_id_start, num_rows);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatPerm<FPP, Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
	{
		return toMatSparse().get_cols(col_ids, num_cols);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatPerm<FPP, Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
	{
		return toMatSparse().get_rows(row_ids, num_rows);
	}

	template<typename FPP>
	std::list<std::pair<int,int>> MatPerm<FPP, Cpu>::nonzeros_indices(const double& tol/*=0*/) const
	{
		auto s = this->getNbRow();
		std::list<std::pair<int,int>> indices;
		const FPP *d_ptr;
		const std::vector<int>* perm_ids_ptr;
		if(is_transp)
		{
			d_ptr = DT.diagonal().data();
			perm_ids_ptr = &perm_ids_T;
		}
		else
		{
			d_ptr = D.diagonal().data();
			perm_ids_ptr = &perm_ids;
		}
		for(int i=0; i < s; i++)
		{
			if(std::abs(d_ptr[i]) > tol)
				indices.push_back(std::make_pair(i, perm_ids_ptr->at(i))); // main diag coeff
		}
		return indices;
	}

	template<typename FPP>
	void MatPerm<FPP, Cpu>::setZeros()
	{
		throw std::runtime_error("setZeros is not available on a MatPerm matrix.");
	}

	template<typename FPP>
	bool MatPerm<FPP, Cpu>::containsNaN() const
	{
		return D.diagonal().hasNaN();
	}

	template<typename FPP>
	const FPP& MatPerm<FPP, Cpu>::operator()(faust_unsigned_int i, faust_unsigned_int j) const
	{
		auto s = this->getNbRow();
		const FPP *d_ptr;
		d_ptr = D.diagonal().data();
		if(i > s || j > s)
			throw std::runtime_error("MatPerm::operator(int i, int j) error: out of bounds coordinates");
		if(j == perm_ids[i])
			return d_ptr[i];
		else
			return FPP(0);
	}


	template<typename FPP>
		MatSparse<FPP, Cpu> MatPerm<FPP, Cpu>::toMatSparse() const
		{
			auto s = this->getNbRow();
			std::vector<Eigen::Triplet<FPP> > tripletList;
			const FPP *d_ptr;
			d_ptr = D.diagonal().data();
			for(int i=0; i < s; i++)
				tripletList.push_back(Eigen::Triplet<FPP>(i, perm_ids[i], d_ptr[i]));
			auto sp = MatSparse<FPP, Cpu>(tripletList, s, s);
			sp.conservativeResize(s, s);
			return sp;
		}


	template<typename FPP>
	bool MatPerm<FPP, Cpu>::isPerm(const MatSparse<FPP, Cpu> &S, bool verify_ones/*=true*/)
	{
		// verify the matrix is square
		if(S.getNbRow() != S.getNbCol()) return false;
		// verify the nnz
		if(S.getNonZeros() != S.getNbRow()) return false;
		// verify only one nz is set per row
		auto rowptr = S.getRowPtr();
		for(int i=0;i < S.getNbRow(); i++)
		{
			if(rowptr[i+1] - rowptr[i] > 1)
				return false;
		}
		// verify all columns are covered
		std::vector<int> cols(S.getNbCol());
		std::iota(cols.begin(), cols.end(), 0);
		auto colind = S.getColInd();
		for(int i=0;i<S.getNonZeros();i++)
		{
			auto it = std::find(cols.begin(), cols.end(), colind[i]);
			if(it == cols.end()) return false; // a column is zero
			cols.erase(it);
		}
		if(cols.size() != 0)
			// not all columns are covered by the nz
			return false;
		if(verify_ones)
		{
			auto valptr = S.getValuePtr();
			for(int i=0;i<S.getNonZeros();i++)
			{
				if(valptr[i] != FPP(1.0))
					return false;
			}
		}
		return true;
	}



	template<typename FPP>
		MatDense<FPP, Cpu> MatPerm<FPP, Cpu>::to_dense() const
		{
			MatDense<FPP, Cpu> dense(this->getNbCol(), this->getNbCol());
			dense.setOnes();
			multiply(dense, 'N');
			return dense;
		}
}
