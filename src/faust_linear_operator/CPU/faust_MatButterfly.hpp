#include "faust_conj.h"
namespace Faust
{

	//TODO: default ctor
	//TODO: dtor?

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
			std::copy(seq.begin()+i+d_offset, seq.begin()+i+2*d_offset, subdiag_ids.begin()+i);
			std::copy(seq.begin()+i, seq.begin()+i+d_offset, subdiag_ids.begin()+i+d_offset);
		}
#ifdef USE_PYTHONIC
		subdiag_ids_ptr = new long[size];
		std::copy(subdiag_ids.begin(), subdiag_ids.end(),subdiag_ids_ptr);
#endif
		this->level = level;
		this->is_transp = false;
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
		std::copy(src.subdiag_ids_ptr, src.subdiag_ids_ptr + src.getNbRow(), subdiag_ids_ptr);
#endif
		level = src.level;
		is_transp = src.is_transp;
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
		std::cout << "subdiag_ids: ";
		for(int i=0;i < subdiag_ids.size();i++)
			std::cout << subdiag_ids[i] << " ";
		std::cout << std::endl;
	}

	template<typename FPP>
		void MatButterfly<FPP, Cpu>::init_transpose()
		{
			//TODO: simplify in case of symmetric matrix (it happens for the FFT)
			if(D2T.size() == 0)
			{
				auto size = D2.rows();
				FPP *d2_ptr, *d2t_ptr;
				d2_ptr = D2.diagonal().data();
				D2T.resize(size);
				d2t_ptr = D2T.diagonal().data();

				auto d_offset = size >> (level+1);
				// D1 doesn't change
				// swap every pair of D2 contiguous blocks to form D2T
				for(int i = 0;i < size; i += d_offset * 2)
				{
					// swap two next blocks of size d_offset into d2t_ptr
					std::copy(d2_ptr + i, d2_ptr + i + d_offset, d2t_ptr + i + d_offset);
					std::copy(d2_ptr + i + d_offset, d2_ptr + i + 2 * d_offset, d2t_ptr + i);
				}
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
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), false);
		return z;
	}


	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(Vect<FPP,Cpu> & x, char opThis) const
	{
		Vect<FPP, Cpu> z(x.size());
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(x.getData(), z.getData(), x.size(), opThis != 'N');
		x = z;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
	{
		MatDense<FPP, Cpu> Y(this->getNbRow(), M.getNbCol());
		const_cast<MatButterfly<FPP, Cpu>*>(this)->multiply(M.getData(), M.getNbCol(), Y.getData(), Y.getNbRow(), 'N' != opThis, 'H' == opThis || 'C' == opThis);
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
	void MatButterfly<FPP, Cpu>::multiply(const FPP* x, FPP* y, size_t size, bool transpose/* = false*/, bool conjugate/* = false*/)
	{
		DiagMat &D2 = this->D2;
		const FPP *d1_ptr, *d2_ptr;
		d1_ptr = D1.diagonal().data();
		auto do_transp = transpose ^ is_transp;
		if(do_transp)
		{
			init_transpose(); // no cost if already initialized
			d2_ptr = D2T.diagonal().data();
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
			if(conjugate)
				y[i] = Faust::conj(d1_ptr[i]) * x[i] + Faust::conj(d2_ptr[i]) * x[subdiag_ids[i]];
			else
				y[i] = d1_ptr[i] * x[i] + d2_ptr[i] * x[subdiag_ids[i]];
#else
		// this is slower
		VecMap x_vec(const_cast<FPP*>(x), size); // const_cast is harmless
		VecMap y_vec(y, size); // const_cast is harmless
		if(do_transp)
			if(conjugate)
				y_vec = diag_conj(D1) * x_vec + diag_conj(D2T) * x_vec(subdiag_ids);
			else
				y_vec = D1 * x_vec + D2T * x_vec(subdiag_ids);
		else
			if(conjugate)
				y_vec = diag_conj(D1) * x_vec + diag_conj(D2) * x_vec(subdiag_ids);
			else
				y_vec = D1 * x_vec + D2 * x_vec(subdiag_ids);
#endif
	}

	template<typename FPP>
		void  MatButterfly<FPP, Cpu>::multiply(const FPP* X, int X_ncols, FPP* Y, size_t Y_nrows, bool transpose/* = false*/, bool conjugate /* = false*/)
		{
			using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
			MatMap X_mat(const_cast<FPP*>(X) /* harmless, no modification*/, Y_nrows, X_ncols);
			MatMap Y_mat(Y, Y_nrows, X_ncols);
			// --------- TODO refactor this block with overload of multiply
			const FPP *d1_ptr, *d2_ptr;
			d1_ptr = D1.diagonal().data();
			auto do_transp = transpose ^ is_transp;
			if(do_transp)
			{
				init_transpose(); // no cost if already initialized
				d2_ptr = D2T.diagonal().data();
			}
			else
				d2_ptr = D2.diagonal().data();
			// ---------------
#if defined(BUTTERFLY_MUL_MAT_OMP_LOOP) || ! EIGEN_VERSION_AT_LEAST(3, 4, 0)
			// this is slower
			#pragma omp parallel for
			for(int i=0;i < Y_nrows; i++)
				if(conjugate)
					Y_mat.row(i) = Faust::conj(d1_ptr[i]) * X_mat.row(i) + Faust::conj(d2_ptr[i]) * X_mat.row(subdiag_ids[i]);
				else
					Y_mat.row(i) = d1_ptr[i] * X_mat.row(i) + d2_ptr[i] * X_mat.row(subdiag_ids[i]);
#else
			// TODO: refactor the exp with a macro diag_prod
			if(do_transp)
				if(conjugate)
					Y_mat = diag_conj(D1) * X_mat + diag_conj(D2T) * X_mat(subdiag_ids, Eigen::placeholders::all);
				else
					Y_mat = D1 * X_mat + D2T * X_mat(subdiag_ids, Eigen::placeholders::all);
			else
				if(conjugate)
					Y_mat = diag_conj(D1) * X_mat + diag_conj(D2) * X_mat(subdiag_ids, Eigen::placeholders::all);
				else
					Y_mat = D1 * X_mat + D2 * X_mat(subdiag_ids, Eigen::placeholders::all);
#endif
		}


	template<typename FPP>
		void MatButterfly<FPP, Cpu>::faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const
		{
			std::runtime_error op_except("Invalid operation letter in MatButterfly::faust_gemm; it must be 'N', 'T' or 'H'.");
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
	void MatButterfly<FPP, Cpu>::transpose()
	{
		init_transpose(); // free cost if already called once
		is_transp = ! is_transp;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::conjugate(const bool eval)
	{
		auto size = getNbRow();
		VecMap d1_vec(const_cast<FPP*>(D1.diagonal().data()), size); // const_cast is harmless
		VecMap d2_vec(const_cast<FPP*>(D2.diagonal().data()), size); // const_cast is harmless
		d1_vec = d1_vec.conjugate();
		d2_vec = d2_vec.conjugate();
		//TODO: D2T
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::adjoint()
	{
		transpose();
		conjugate();
	}

	template<typename FPP>
	faust_unsigned_int MatButterfly<FPP, Cpu>::getNonZeros()const
	{
		return D1.diagonal().nonZeros() + D2.diagonal().nonZeros();
	}

	template<typename FPP>
	size_t MatButterfly<FPP, Cpu>::getNBytes() const
	{
		return (D1.rows() + D2.rows() + (D2T.size() != 0?D2T.rows():0)) * sizeof(FPP) + subdiag_ids.size() * sizeof(int);
	}

	template<typename FPP>
	MatType MatButterfly<FPP, Cpu>::getType() const
	{
		return Butterfly;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::operator*=(const FPP alpha)
	{
		D1 = D1 * alpha;
		D2 = D2 * alpha;
		if(D2T.size() != 0)
			D2T = D2T * alpha;
	}

	template<typename FPP>
	matvar_t* MatButterfly<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate, const char *var_name) const
	{
		throw std::runtime_error("Saving MatButterfly to a .mat file is not supported.");
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::normL1(const bool transpose/*=false*/) const
	{
		return toMatSparse().normL1(transpose);
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::norm() const
	{
		auto s = (D1.diagonal().array() * D1.diagonal().conjugate().array() + D2.diagonal().array() * D2.diagonal().conjugate().array()).sum();
		return Faust::fabs(std::sqrt(s));
	}

	template<typename FPP>
	Real<FPP> MatButterfly<FPP, Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose) const
	{
		return toMatSparse().normL1(col_id, transpose);
	}

	template<typename FPP>
	Vect<FPP,Cpu> MatButterfly<FPP, Cpu>::get_col(faust_unsigned_int id) const
	{
		return toMatSparse().get_col(id);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
	{
		return toMatSparse().get_cols(col_id_start, num_cols);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
	{
		return toMatSparse().get_rows(row_id_start, num_rows);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
	{
		return toMatSparse().get_cols(col_ids, num_cols);
	}

	template<typename FPP>
	MatGeneric<FPP,Cpu>* MatButterfly<FPP, Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
	{
		return toMatSparse().get_rows(row_ids, num_rows);
	}

	template<typename FPP>
	std::list<std::pair<int,int>> MatButterfly<FPP, Cpu>::nonzeros_indices(const double& tol/*=0*/) const
	{
		//TODO: transpose case
		auto s = this->getNbRow();
		auto k = s >> (level + 1); // D2 diagonal offset
		std::list<std::pair<int,int>> indices;
		const FPP *d1_ptr, *d2_ptr;
		d1_ptr = D1.diagonal().data();
		d2_ptr = D2.diagonal().data();
		for(int i=0; i < s; i++)
		{
			if(std::abs(d1_ptr[i]) > tol)
				indices.push_back(std::make_pair(i, i)); // main diag coeff
			if(std::abs(d2_ptr[i]) > tol)
				if((i / k) & 1)
					indices.push_back(std::make_pair(i, i - k)); // below diag coeff
				else
					indices.push_back(std::make_pair(i, i + k)); // above diag coeff
		}
		return indices;
	}

	template<typename FPP>
	void MatButterfly<FPP, Cpu>::setZeros()
	{
		throw std::runtime_error("setZeros is not available on a MatButterfly matrix.");
	}

	template<typename FPP>
	bool MatButterfly<FPP, Cpu>::containsNaN() const
	{
		return D1.diagonal().hasNaN() || D2.diagonal().hasNaN();
	}

	template<typename FPP>
	const FPP& MatButterfly<FPP, Cpu>::operator()(faust_unsigned_int i, faust_unsigned_int j) const
	{
		auto s = this->getNbRow();
		if(i > s || j > s)
			throw std::runtime_error("MatButterfly::operator(int i, int j) error: out of bounds coordinates");
		auto k = s >> (level + 1); // D2 diagonal offset
		const FPP *d1_ptr, *d2_ptr;
		d1_ptr = D1.diagonal().data();
		d2_ptr = D2.diagonal().data();
		if(i == j)
			return d1_ptr[i];
		if((i / k) & 1)
			if(j == i - k)
				return d2_ptr[i];
			else
				return FPP(0);
		if(j == i + k)
			return d2_ptr[i];
		else
			return FPP(0);
	}


	template<typename FPP>
		MatSparse<FPP, Cpu> MatButterfly<FPP, Cpu>::toMatSparse() const
		{
			auto s = this->getNbRow();
			auto k = s >> (level + 1); // D2 diagonal offset
			std::vector<Eigen::Triplet<FPP> > tripletList;
			const FPP *d1_ptr, *d2_ptr;
			d1_ptr = D1.diagonal().data();
			d2_ptr = D2.diagonal().data();
			for(int i=0; i < s; i++)
			{
				tripletList.push_back(Eigen::Triplet<FPP>(i, i, d1_ptr[i]));
				if((i / k) & 1)
				{
					// d2 coeff is the first elt of row i
					tripletList.push_back(Eigen::Triplet<FPP>(i, i - k, d2_ptr[i]));
				}
				else
				{
					// d2 coeff is the last elt of row i
					tripletList.push_back(Eigen::Triplet<FPP>(i, i + k, d2_ptr[i]));
				}
			}
			auto sp = MatSparse<FPP, Cpu>(tripletList, s, s);
			sp.conservativeResize(s, s);
			return sp;
		}

}
