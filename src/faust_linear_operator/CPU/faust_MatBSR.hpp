#include <stdexcept>
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <numeric>


namespace Faust
{
	template<typename FPP>
		MatBSR<FPP, Cpu>::MatBSR(faust_unsigned_int nrows, faust_unsigned_int ncols, faust_unsigned_int bnrows, faust_unsigned_int bncols, faust_unsigned_int nblocks, const FPP* data, const int *block_rowptr, const int *block_colinds) : MatGeneric<FPP, Cpu>(nrows, ncols), bmat(nrows, ncols, bnrows, bncols, nblocks, data, block_rowptr, block_colinds)
	{
	}

	template <typename FPP>
		MatBSR<FPP,Cpu>::MatBSR(BSRMat<FPP>& bmat)
		{
			this->bmat = bmat;
			this->dim1 = bmat.m;
			this->dim2 = bmat.n;
		}

	template <typename FPP>
		MatGeneric<FPP,Cpu>* MatBSR<FPP,Cpu>::Clone(const bool isOptimize/*=false*/) const
		{
			//TODO: optimize to a MatDense or MatSparse by comparing the time of matrix-vector product?
			MatBSR<FPP,Cpu> *clone = new MatBSR<FPP,Cpu>();
			BSRMat<FPP> clone_bmat(bmat);
			clone->bmat = clone_bmat;
			// clone->bmat.print_bufs();
			clone->dim1 = clone_bmat.m;
			clone->dim2 = clone_bmat.n;
			return clone;
		}


	template <typename FPP>
		void MatBSR<FPP,Cpu>::multiply(Faust::Vect<FPP,Cpu> & vec, char opThis) const
		{
			switch(opThis)
			{
				case 'N': {
					vec.vec = bmat.mul(vec.vec);
						  }
				break;
				case 'T': {
					auto bmat_t = bmat;
					bmat_t.transpose(/*inplace*/ true);
					vec.vec = bmat_t.mul(vec.vec);
						  }
				break;
				case 'H': {
					auto bmat_tc = bmat;
					bmat_tc.transpose(/*inplace*/ true);
					bmat_tc.conjugate(/*inplace*/ true);
					vec.vec = bmat_tc.mul(vec.vec);
						  }
				break;
				default:
				{
					throw std::runtime_error("Unknown op type.");
				}
			}
		}

	template <typename FPP>
	Vect<FPP,Cpu> MatBSR<FPP,Cpu>::multiply(const Vect<FPP,Cpu> &v) const
	{
		Vect<FPP, Cpu> vec;
		vec.resize(this->dim1);
		vec.vec = bmat.mul(v.vec);
		return vec;
	}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::multiply(Faust::MatDense<FPP,Cpu> & M, char opThis) const
		{
			switch(opThis)
			{
				case 'N':
					{
						M.mat = bmat.mul(M.mat);
						break;
					}
				case 'T':
					{
						auto bmat_t = bmat;
						bmat_t.transpose(/*inplace*/ true);
						M.mat = bmat_t.mul(M.mat);
						break;
					}
				case 'H':
					{
						auto bmat_tc = bmat;
						bmat_tc.transpose(/*inplace*/ true);
						bmat_tc.conjugate(/*inplace*/ true);
						M.mat = bmat_tc.mul(M.mat);
						break;
					}
				default:
					{
						throw std::runtime_error("Unknown op type.");
					}
			}
			M.update_dims();
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::multiply(Faust::MatSparse<FPP,Cpu> & M, char opThis) const
		{
			switch(opThis)
			{
				case 'N':
					{
						M.mat = bmat.mul(M.mat).sparseView();
						break;
					}
				case 'T':
					{
						auto bmat_t = bmat;
						bmat_t.transpose(/*inplace*/ true);
						M.mat = bmat_t.mul(M.mat).sparseView();
						break;
					}
				case 'H':
					{
						auto bmat_tc = bmat;
						bmat_tc.transpose(/*inplace*/ true);
						bmat_tc.conjugate(/*inplace*/ true);
						M.mat = bmat_tc.mul(M.mat).sparseView();
						break;
					}
				default:
					{
						throw std::runtime_error("Unknown op type.");
					}
			}
			M.update_dim();
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::multiplyRight(Faust::MatSparse<FPP, Cpu> const& M)
		{
			// normally at BSRMat times a MatSparse gives a DenseMat
			// so the only way to return a MatBSR is to construct a MatBSR (in fact a BSRMat) composed of one single block which is the whole matrix
			// However this matrix has most likely no interest to be a MatBSR and should be converted to a MatDense
			auto dmat = bmat.mul(M.mat);
			BSRMat<FPP> bmat;
			size_t data_size = dmat.rows()*dmat.cols();
			bmat.data = new FPP[data_size];
			memcpy(bmat.data, dmat.data(), sizeof(FPP)*data_size);
			bmat.bcolinds = new int[1];
			bmat.bcolinds[0] = 0;
			bmat.browptr = new int[2];
			bmat.browptr[0] = 0;
			bmat.browptr[1] = 1;
			this->dim1 = bmat.m = dmat.rows();
			this->dim2 = bmat.n = dmat.cols();
			bmat.bnnz = 1;
			bmat.bm = dmat.rows();
			bmat.bn = dmat.cols();
			bmat.b_per_rowdim = bmat.b_per_coldim = 1;
			this->bmat = bmat; // operator= frees bmat buffers before override
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const // from LinearOperator
		{
			MatBSR<FPP, Cpu> this_copy;
			MatDense<FPP, Cpu> B_copy(B);
			if(typeB == 'T')
				B_copy.transpose();
			else if(typeB == 'H')
				B_copy.adjoint();
			if(B.getNBytes() > this->getNBytes())
			{
				// B is heavier than this
				// do the scalar mul on this
				this_copy = *this;
				this_copy.bmat.mul(alpha);
				this_copy.multiply(B_copy, typeA);
			}
			else
			{
				B_copy *= alpha;
				multiply(B_copy, typeA);
			}
			C *= beta;
			C.add(B_copy);
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::transpose()
		{
			bmat.transpose(/*inplace*/ true);
			this->dim1 = bmat.m;
			this->dim2 = bmat.n;
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::conjugate(const bool eval/*=true*/)
		{
			bmat.conjugate(/*inplace*/ true);
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::adjoint()
		{
			bmat.adjoint(/*inplace*/ true);
			this->dim1 = bmat.m;
			this->dim2 = bmat.n;
		}

	template <typename FPP>
		faust_unsigned_int MatBSR<FPP,Cpu>::getNonZeros()const
		{
			return bmat.nnz();
		}

	template <typename FPP>
		size_t MatBSR<FPP,Cpu>::getNBytes()const
		{
			return bmat.nbytes();
		}

	template <typename FPP>
		size_t MatBSR<FPP,Cpu>::getNbBlockRow()const
		{
			return bmat.bm;
		}

	template <typename FPP>
		size_t MatBSR<FPP,Cpu>::getNbBlockCol()const
		{
			return bmat.bn;
		}

	template <typename FPP>
		size_t MatBSR<FPP,Cpu>::getNbBlocksPerDim(int dim_id) const
		{
			switch(dim_id)
			{
				case 0:
					return bmat.b_per_rowdim;
				case 1:
					return bmat.b_per_coldim;
				default:
					throw std::runtime_error("MatBSR::getNbBlocksPerDim invalid dimension id (must be 0 or 1).");
			}
		}

	template <typename FPP>
		size_t MatBSR<FPP,Cpu>::getNBlocks() const
		{
			return bmat.bnnz;
		}

	template <typename FPP>
		MatType MatBSR<FPP,Cpu>::getType() const
		{
			return BSR;
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::operator*=(const FPP alpha)
		{
			bmat.mul(alpha);
		}

	template <typename FPP>
		std::string MatBSR<FPP,Cpu>::to_string(MatType type, const bool transpose/*=false*/, const bool displaying_small_mat_elts/*=false*/) const
		{
			std::ostringstream  str;
			str << MatGeneric<FPP,Cpu>::to_string(BSR, transpose);
			if(bmat.data == nullptr)
				str <<"zeros matrix flag" <<std::endl;
			else
			{
				auto dmat = bmat.to_dense();
				if (displaying_small_mat_elts && this->dim1*this->dim2 < 1000)
				{
					for (int i=0 ; i<this->dim1 ; i++)
					{
						for(int j=0 ; j<this->dim2 ; j++)
							str << dmat(i,j) << " " ;
						str << std::endl;
					}
				}
			}
			return str.str();
		}

	template <typename FPP>
		std::string MatBSR<FPP,Cpu>::to_string(const bool transpose/*=false*/, const bool displaying_small_mat_elts/*=false*/) const
		{
			return MatGeneric<FPP,Cpu>::to_string(this->getNbRow(), this->getNbCol(), transpose, bmat.density(), getNonZeros(), this->is_identity, BSR);
		}

	template <typename FPP>
		matvar_t* MatBSR<FPP,Cpu>::toMatIOVar(bool transpose, bool conjugate) const
		{
			//			MatSparse<FPP, Cpu> smat(this->dim1, this->dim2);
			//			smat.mat = bmat.to_sparse();
			//			return smat.toMatIOVar(transpose, conjugate);
			int opt = typeid(bmat.data[0])==typeid(std::complex<Real<FPP>>(1.0,1.0))?MAT_F_COMPLEX:0;
			size_t bsr_name_dims[2] = {1, 3};
			const char* bsr_name = "bsr";
			size_t params_dims[2] = {1, 3}; // nrows, ncols, bnnz
			int params[3] = {(int) this->getNbRow(), (int) this->getNbCol(), (int) this->getNBlocks()};
			size_t nzb_dims[2] = {1, getNBlocks()}; // bcolinds size
			size_t browptr_dims[2] = {1, getNbBlocksPerDim(0)+1};
			size_t bcolinds_dims[2] = {1, getNBlocks()};
			size_t bdata_dims[2] = {1, getNBlocks()*getNbBlockRow()*getNbBlockCol()};
			size_t cell_dims[2] = { 1, 5};
			matvar_t *matvars[6], *bsr_cell = nullptr;
			mat_complex_split_t z = {nullptr,nullptr};
			matio_types bdata_matio_type;
			matio_classes bdata_matio_class;
			using DenseMatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
			if(std::is_same<FPP, float>::value)
			{
				bdata_matio_type = MAT_T_SINGLE;
				bdata_matio_class = MAT_C_SINGLE;
			}
			else
			{
				bdata_matio_type = MAT_T_DOUBLE;
				bdata_matio_class = MAT_C_DOUBLE;
			}
			matvars[5] = nullptr;
			matvars[0] = Mat_VarCreate("bsr_name", MAT_C_CHAR, MAT_T_UINT8, 2, bsr_name_dims, (void*) bsr_name, 0); // copy needed because bsr_name is in the stack
//			auto mat_bsr_name = Mat_CreateVer("test_bsr_name.mat", NULL, MAT_FT_MAT5);
//			Mat_VarWrite(mat_bsr_name, matvars[0], MAT_COMPRESSION_NONE);
//			Mat_Close(mat_bsr_name);
			matvars[1] = Mat_VarCreate("params", MAT_C_INT32, MAT_T_INT32, 2, params_dims, params, 0); // copy needed because params is in the stack
			matvars[2] = Mat_VarCreate("bcolinds", MAT_C_INT32, MAT_T_INT32, 2, bcolinds_dims, bmat.bcolinds, MAT_F_DONT_COPY_DATA);
			matvars[3] = Mat_VarCreate("browptr", MAT_C_INT32, MAT_T_INT32, 2, browptr_dims, bmat.browptr, MAT_F_DONT_COPY_DATA);
//			auto mat_browptr = Mat_CreateVer("test_bsr_browptr.mat", NULL, MAT_FT_MAT5);
//			Mat_VarWrite(mat_browptr, matvars[3], MAT_COMPRESSION_NONE);
//			Mat_Close(mat_browptr);
//			auto mat_bcolinds = Mat_CreateVer("test_bsr_bcolinds.mat", NULL, MAT_FT_MAT5);
//			Mat_VarWrite(mat_bcolinds, matvars[2], MAT_COMPRESSION_NONE);
//			Mat_Close(mat_bcolinds);
			DenseMatMap _bmat(bmat.data, getNbBlockRow(), getNBlocks()*getNbBlockCol());
			if(opt)
			{
				// complex matrix
				DenseMat<Real<FPP>> dst_re(getNbBlockRow(), getNBlocks()*getNbBlockCol());
				DenseMat<Real<FPP>> dst_im(getNbBlockRow(), getNBlocks()*getNbBlockCol());
//				DenseMat<FPP> tmp(getNbBlockRow(), getNBlocks()*getNbBlockCol());
//				tmp = _bmat;
				if(transpose)
				{
					dst_re = _bmat.real().template cast<Real<FPP>>().transpose();
//					dst_im = _bmat.imag().template cast<Real<FPP>>().transpose();
					for(int i=0;i<_bmat.rows()*_bmat.cols();i++)
						dst_im.data()[i] = std::imag(_bmat.data()[i]);
					dst_im.transposeInPlace();
				}
				else
				{
					dst_re = _bmat.real().template cast<Real<FPP>>();
//					dst_im = _bmat.imag().template cast<Real<FPP>>();
					for(int i=0;i<_bmat.rows()*_bmat.cols();i++)
						dst_im.data()[i] = std::imag(_bmat.data()[i]);

				}
				if(conjugate)
					dst_im *= Real<FPP>(-1);
				z.Re = dst_re.data();
				z.Im = dst_im.data();
				matvars[4] = Mat_VarCreate(nullptr, bdata_matio_class, bdata_matio_type, 2, bdata_dims, &z /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);

			}
			else
			{
				// real matrix // conjugate is ignored
				if(transpose)
				{
					DenseMat<FPP> mat(getNbBlockRow(), getNBlocks()*getNbBlockCol());
					mat = _bmat.transpose();
					matvars[4] = Mat_VarCreate("bdata", bdata_matio_class, bdata_matio_type, 2, bdata_dims, mat.data(), MAT_F_DONT_COPY_DATA);
				}
				else
					matvars[4] = Mat_VarCreate("bdata", bdata_matio_class, bdata_matio_type, 2, bdata_dims, bmat.data, MAT_F_DONT_COPY_DATA);
			}
//			auto mat = Mat_CreateVer("test_bsr_bdata.mat", NULL, MAT_FT_MAT5);
//			Mat_VarWrite(mat, matvars[4], MAT_COMPRESSION_NONE);
//			Mat_Close(mat);
			bsr_cell = Mat_VarCreate("bsr_cell", MAT_C_CELL, MAT_T_CELL, 2, cell_dims, matvars, 0);
			return bsr_cell;
		}

	template <typename FPP>
		Real<FPP> MatBSR<FPP,Cpu>::normL1(const bool transpose/*=false*/) const
		{
			return bmat.normL1();
		}

	template <typename FPP>
		Real<FPP> MatBSR<FPP,Cpu>::norm() const
		{
			return bmat.norm();
		}

	template <typename FPP>
		Real<FPP> MatBSR<FPP,Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose) const
		{
			//TODO:
			throw std::runtime_error("Not supported yet.");
		}

	template <typename FPP>
		Vect<FPP,Cpu> MatBSR<FPP,Cpu>::get_col(faust_unsigned_int id) const
		{
			Vect<FPP, Cpu> v(this->dim1);
			v.vec = bmat.get_col(id);
			return v;
		}

	template <typename FPP>
		MatSparse<FPP,Cpu>* MatBSR<FPP,Cpu>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
		{
			auto smat = new MatSparse<FPP, Cpu>(this->dim1, num_cols);
			smat->mat = bmat.get_cols(col_id_start, num_cols);
			smat->update_dim();
			return smat;
		}

	template <typename FPP>
		MatSparse<FPP,Cpu>* MatBSR<FPP,Cpu>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
		{
			auto smat = new MatSparse<FPP, Cpu>(num_rows, this->dim2);
			smat->mat = bmat.get_rows(row_id_start, num_rows);
			smat->update_dim();
			return smat;
		}

	template <typename FPP>
		MatSparse<FPP,Cpu>* MatBSR<FPP,Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
		{
			auto smat = new MatSparse<FPP, Cpu>(this->dim1, num_cols);
			smat->mat = bmat.get_cols(col_ids, num_cols);
			smat->update_dim();
			return smat;
		}

	template <typename FPP>
		MatSparse<FPP,Cpu>* MatBSR<FPP,Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
		{
			auto smat = new MatSparse<FPP, Cpu>(num_rows, this->dim2);
			smat->mat = bmat.get_rows(row_ids, num_rows);
			smat->update_dim();
			return smat;
		}

	template <typename FPP>
		std::list<std::pair<int,int>> MatBSR<FPP,Cpu>::nonzeros_indices() const
		{
			return bmat.nonzeros_indices();
		}

	template <typename FPP>
		void MatBSR<FPP,Cpu>::setZeros()
		{
			//TODO: verify how the operations behave after setZeros
			bmat.free_bufs();
			bmat.bnnz = 0;
		}

	template <typename FPP>
		bool MatBSR<FPP,Cpu>::containsNaN() const
		{
			return bmat.contains_nan();
		}

	template <typename FPP>
		const FPP& MatBSR<FPP,Cpu>::operator()(faust_unsigned_int i, faust_unsigned_int j)const
		{
			return bmat(i,j);
		}

	template <typename FPP>
		MatDense<FPP, Cpu> MatBSR<FPP,Cpu>::to_dense()const
		{
			MatDense<FPP, Cpu> dmat;
			dmat.mat = bmat.to_dense();
			dmat.update_dims();
			return dmat;
		}

	template <typename FPP>
		MatSparse<FPP, Cpu> MatBSR<FPP,Cpu>::to_sparse()const
		{
			MatSparse<FPP, Cpu> smat(this->dim1, this->dim2);
			smat.mat = bmat.to_sparse();
			smat.update_dim();
			return smat;
		}

	template <typename FPP>
		MatBSR<FPP,Cpu>::~MatBSR()
		{
		}

	template <typename FPP>
		MatBSR<FPP,Cpu>* MatBSR<FPP,Cpu>::randMat(int m, int n, int bm, int bn, int bnnz)
		{
			auto bmat = BSRMat<FPP>::rand(m,n,bm,bn,bnnz);
			auto rbmat = new MatBSR<FPP,Cpu>(bmat);
			return rbmat;
		}

}


/******************** Low-level BSRMat definitions ***********************/
template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>::BSRMat(const BSRMat<T, BlockStorageOrder>& src_bmat) : BSRMat<T, BlockStorageOrder>()
{
	*this = src_bmat;
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>::BSRMat(unsigned long int nrows, unsigned long int ncols, unsigned long int bnrows, unsigned long int bncols, unsigned long int nblocks, const T* data, const int *block_rowptr, const int *block_colinds) : BSRMat<T, BlockStorageOrder>()
{
		// verify bm and bn evenly divide m and n
		if(nrows%bnrows)
			throw std::runtime_error("BSRMat error: bnrows must evenly divide nrows.");
		if(ncols%bncols)
			throw std::runtime_error("BSRMat error: bncols must evenly divide ncols.");
		// enforce bnnz is less or equel to m*n/bm/bn
		nblocks = std::min(nblocks, nrows*ncols/bnrows/bncols);
		this->m = nrows;
		this->n = ncols;
		this->bm = bnrows;
		this->bn = bncols;
		this->bnnz = nblocks;
		this->b_per_rowdim = this->m/this->bm;
		this->b_per_coldim = this->n/this->bn;
		int data_size = bnnz*bm*bn;
		// init data
		this->data = new T[data_size];
		memcpy(this->data, data, sizeof(T)*nblocks*bm*bn);
		this->browptr = new int[this->b_per_rowdim+1];
		memcpy(this->browptr, block_rowptr, sizeof(int)*(this->b_per_rowdim+1));
		this->bcolinds = new int[this->bnnz];
		memcpy(this->bcolinds, block_colinds, sizeof(int)*this->bnnz);
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>::BSRMat(BSRMat<T, BlockStorageOrder>&& src_bmat) : BSRMat<T, BlockStorageOrder>()
{
	*this = std::move(src_bmat);
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>& BSRMat<T, BlockStorageOrder>::operator=(BSRMat<T, BlockStorageOrder>&& src_bmat)
{
	bnnz = src_bmat.bnnz;
	m = src_bmat.m;
	n = src_bmat.n;
	bm = src_bmat.bm;
	bn = src_bmat.bn;
	b_per_rowdim = src_bmat.b_per_rowdim;
	b_per_coldim = src_bmat.b_per_coldim;
	data = src_bmat.data;
	browptr = src_bmat.browptr;
	bcolinds = src_bmat.bcolinds;
	src_bmat.data = nullptr;
	src_bmat.browptr = src_bmat.bcolinds = nullptr;
	return *this;
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>& BSRMat<T, BlockStorageOrder>::operator=(const BSRMat<T, BlockStorageOrder>& src_bmat)
{
	bnnz = src_bmat.bnnz;
	m = src_bmat.m;
	n = src_bmat.n;
	bm = src_bmat.bm;
	bn = src_bmat.bn;
	b_per_rowdim = src_bmat.b_per_rowdim;
	b_per_coldim = src_bmat.b_per_coldim;
	free_bufs();
	data = new T[bnnz*bm*bn];
	bcolinds = new int[bnnz];
	browptr = new int[b_per_coldim+1];
	memcpy(data, src_bmat.data, sizeof(T)*bnnz*bm*bn);
	memcpy(bcolinds, src_bmat.bcolinds, sizeof(int)*bnnz);
	memcpy(browptr, src_bmat.browptr, sizeof(int)*(b_per_coldim+1));
	return *this;
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder> BSRMat<T, BlockStorageOrder>::rand(int m, int n, int bm, int bn, int bnnz)
{
	// verify bm and bn evenly divide m and n
	if(m%bm)
		throw std::runtime_error("BSRMat::rand error: bm must evenly divide m.");
	if(n%bn)
		throw std::runtime_error("BSRMat::rand error: bn must evenly divide n.");
	// enforce bnnz is less or equel to m*n/bm/bn
	bnnz = std::min(bnnz, m*n/bm/bn);
	// instantiate a BSRMat
	BSRMat bmat;
	bmat.m = m;
	bmat.n = n;
	bmat.bm = bm;
	bmat.bn = bn;
	bmat.bnnz = bnnz;
	bmat.b_per_rowdim = m/bm;
	bmat.b_per_coldim = n/bn;
	auto nblocks = bmat.b_per_coldim*bmat.b_per_rowdim;
	// choose bnnz (i,j) indices to init bcolinds and browptr
	std::random_device rd;
	std::mt19937 gen(rd());
	// data size ?
	int data_size = bnnz*bm*bn;
	// init data
	bmat.data = new T[data_size];
	std::uniform_real_distribution<Real<T>> dis(0.0, 1.0);
	for(int i=0;i<data_size;i++)
		bmat.data[i] = dis(gen);
	// generate the vector of possible block indices (in 1D first to ease the selection of bnnz random blocks)
	std::vector<int> BI(nblocks);
	std::vector<int> nnzBI(bnnz);
	std::iota(BI.begin(), BI.end(), 0);
	// vector of selected 1D indices of bnnz blocks
	for(int i = 0; i < bnnz; i++)
	{
		std::uniform_int_distribution<> distrib_block_inds(0, BI.size()-1);
		int r = distrib_block_inds(gen);
		int bi = BI[r];
		BI.erase(BI.begin()+r);
		nnzBI[i] = bi;
	}
	std::sort(nnzBI.begin(), nnzBI.end());
	bmat.browptr = new int[bmat.b_per_rowdim+1];
	bmat.bcolinds = new int[bnnz];
	bmat.browptr[0] = 0;
	int browptr_ind = 1;
	int brow_nblocks = 0;
	int brow_ind = 0;
	int bcolinds_ind = 0;
	memset(bmat.browptr, 0, sizeof(int)*(bmat.b_per_rowdim+1));
	for(auto bi: nnzBI)
	{
		while(bi/bmat.b_per_coldim+1 > browptr_ind)
			browptr_ind++;
		// auto rowblock_ind = bi/bn;
		bmat.bcolinds[bcolinds_ind++] = bi%bmat.b_per_coldim;
		bmat.browptr[browptr_ind]++;
	}
	for(int i=1; i < bmat.b_per_rowdim+1; i++)
		bmat.browptr[i] += bmat.browptr[i-1];
	return bmat;
}


template<typename T, int BlockStorageOrder>
void BSRMat<T,BlockStorageOrder>::print_bufs()
{
	int block_size = bm*bn;
	int data_size = bnnz*block_size;
	std::cout << "data (nz block-by-block):" << std::endl;
	for(int i = 0; i < data_size;i++)
	{
		std::cout << data[i] << " ";
		if((i+1)%block_size == 0)
			std::cout << std::endl;
	}
	std::cout << "cumulative num of blocks / blocks per row:" << std::endl;
	for(int i=1;i<m/bm+1;i++)
	{
		std::cout << browptr[i] << " " << browptr[i]-browptr[i-1] << std::endl;
	}
	std::cout << "bcolinds" << std::endl;
	for(int j=0;j<bnnz;j++)
	{
		std::cout << bcolinds[j] << " ";
	}
	std::cout << std::endl;
}

template<typename T, int BlockStorageOrder> DenseMat<T> BSRMat<T, BlockStorageOrder>::to_dense() const
{
	DenseMat<T> dmat = DenseMat<T>::Zero(this->m, this->n);
	iter_block([&dmat, this](int mat_row_id, int mat_col_id, int block_offset)
			{
				dmat.block(mat_row_id, mat_col_id, bm, bn) = DenseMatMap<T>(data+block_offset*bm*bn, bm, bn);
			});
	return dmat;
}

template<typename T, int BlockStorageOrder> SparseMat<T> BSRMat<T, BlockStorageOrder>::to_sparse() const
{
	SparseMat<T> smat(this->m, this->n);
	typedef Eigen::Triplet<T,int> TRI;
	std::vector<TRI> triplets;
	triplets.reserve(bnnz*bm*bn);
	iter_block([&smat, this, &triplets](int mat_row_id, int mat_col_id, int block_offset)
			{
			for(int i=0;i<bm;i++)
				for(int j=0;j<bn;j++)
						triplets.push_back(Eigen::Triplet<T>(mat_row_id+i, mat_col_id+j, (T)data[block_offset*bm*bn+j*bm+i]));
			});
	smat.setFromTriplets(triplets.begin(), triplets.end());
	smat.makeCompressed();
	return smat;
}

template<typename T, int BlockStorageOrder>
Vec<T> BSRMat<T, BlockStorageOrder>::mul(const T* vec_data) const
{
	Vec<T> prod = Vec<T>::Zero(m, 1);
	iter_block([&prod, vec_data, this](int mat_row_id, int mat_col_id, int block_offset)
			{
				prod.block(mat_row_id, 0, bm, 1) += DenseMatMap<T>(data+block_offset*bm*bn, bm, bn) * DenseMatMap<T>(vec_data+mat_col_id, bn, 1);
			});
	return prod;
}

template<typename T, int BlockStorageOrder>
DenseMat<T> BSRMat<T, BlockStorageOrder>::mul(const DenseMat<T>& other) const
{
	DenseMat<T> prod = DenseMat<T>::Zero(m, other.cols());
	iter_block([&prod, &other, this](int mat_row_id, int mat_col_id, int block_offset)
			{
			prod.block(mat_row_id, 0, bm, other.cols()) += DenseMatMap<T>(data+block_offset*bm*bn, bm, bn) * other.block(mat_col_id, 0, bn, other.cols());
			});
	return prod;
}

template<typename T, int BlockStorageOrder>
DenseMat<T> BSRMat<T, BlockStorageOrder>::mul(const SparseMat<T>& other) const
{
	DenseMat<T> prod = DenseMat<T>::Zero(m, other.cols());
	iter_block([&prod, &other, this](int mat_row_id, int mat_col_id, int block_offset)
			{
			prod.block(mat_row_id, 0, bm, other.cols()) += DenseMatMap<T>(data+block_offset*bm*bn, bm, bn) * other.block(mat_col_id, 0, bn, other.cols());
			});
	return prod;
}

template<typename T, int BlockStorageOrder>
DenseMat<T> BSRMat<T, BlockStorageOrder>::mul(const BSRMat<T, BlockStorageOrder>& other) const
{
	auto smat = other.to_sparse();
	return mul(smat);
}

template<typename T, int BlockStorageOrder>
void BSRMat<T, BlockStorageOrder>::mul(const T& scal) const
{
	DenseMatMap<T>(data, bm*bn*bnnz, 1) *= scal;
}

template<typename T, int BlockStorageOrder>
void BSRMat<T, BlockStorageOrder>::iter_block(std::function<void(int /* mat_row_id */, int /* mat_col_id */, int /* block_offset */)> op) const
{
	int mat_row_id, mat_col_id;
	int block_offset = 0;
	int i = 1; // current line of blocks + 1
	int j = browptr[i]; // the number of nonzero blocks at the current line of blocks
	// iterate on each row of blocks and call op for each block
	while(i < this->b_per_rowdim+1 && block_offset < bnnz)
	{
		// find the next nonzero line of blocks
		while(j == 0 && i < this->b_per_rowdim+1)
		{
			i++;
			if(i < this->b_per_rowdim+1)
				j = browptr[i]-browptr[i-1];
		}
		if(i < this->b_per_rowdim+1)
		{
			mat_row_id = (i-1)*bm;
			mat_col_id = bcolinds[block_offset]*bn;
			op(mat_row_id, mat_col_id, block_offset);
			assert(block_offset < bnnz);
			block_offset++;
			j--;
		}
	}
}

template<typename T, int BlockStorageOrder>
Real<T> BSRMat<T, BlockStorageOrder>::norm() const
{
	return DenseMatMap<T>(data, bm*bnnz, bn).norm();
}

template<typename T, int BlockStorageOrder>
Real<T> BSRMat<T, BlockStorageOrder>::normInf() const
{
	Real<T> maxSum = 0;
	Real<T> *sums = new Real<T>[bm];
	memset(sums, 0, sizeof(Real<T>)*bm);
	int cur_row = 0;
	auto _max_sums = [this, &sums, &maxSum](int cur_row)
	{
		for(int i=0; i<bm; i++)
		{
//			std::cout << "sum row i=" << cur_row+i << ": " << sums[i] << std::endl;
			if(sums[i] > maxSum)
				maxSum = sums[i];
		}
	};
	iter_block([&maxSum, &sums, &cur_row, this, _max_sums](int mat_row_id, int mat_col_id, int block_offset)
			{
				if(cur_row < mat_row_id)
				{
					// new block line
					// compute the inf-norm for the previous block lines
					_max_sums(cur_row);
					cur_row = mat_row_id;
					memset(sums, 0, sizeof(Real<T>)*bm);
				}
//				std::cout << "block_offset:" << block_offset << std::endl;
				for(int i=0; i<bm; i++)
				{
//					std::cout << "row: " << cur_row+i << ": ";
					for(int j=0; j<bn; j++)
					{
//						std::cout << data[block_offset*bm*bn+j*bm+i] << " ";
						sums[i]	+= std::abs(data[block_offset*bm*bn+j*bm+i]);
					}
//					std::cout << std::endl;
				}
			});
	_max_sums((b_per_rowdim-1)*bm);
	delete[] sums;
	return maxSum;
}

template<typename T, int BlockStorageOrder>
Real<T> BSRMat<T, BlockStorageOrder>::normL1() const
{
	auto t_this = *this;
	t_this.transpose(/*inplace*/true);
	return t_this.normInf();
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder> BSRMat<T, BlockStorageOrder>::transpose(const bool inplace/*=false*/)
{
	return apply_op('T', inplace);
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder> BSRMat<T, BlockStorageOrder>::conjugate(const bool inplace/*="false"*/)
{
	if(inplace)
	{
		DenseMatMap<T> data_mat(data, bnnz*bm, bn);
		data_mat = DenseMatMap<T>(data, bnnz*bm, bn).conjugate();
		return *this;
	}
	else
	{
		BSRMat<T, BlockStorageOrder> cbmat(*this);
		cbmat.conjugate(/*inplace=*/true);
		return cbmat;
	}
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder> BSRMat<T, BlockStorageOrder>::adjoint(const bool inplace/*="false"*/)
{
	return apply_op('H', inplace);
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder> BSRMat<T, BlockStorageOrder>::apply_op(const char op, const bool inplace/*="false"*/)
{
	if(op != 'C' && op != 'H' && op != 'T' && op != 'N')
		throw std::runtime_error("BSRMat::apply_op: unknown op.");
	if(op == 'N')
		if(inplace)
			return *this;
		else
			return BSRMat<T, BlockStorageOrder>(*this);
	if(op == 'C')
		return conjugate(inplace);
	// op == 'T' or 'H'
	if(inplace)
	{
		// init tbmat buffers and attributes
		auto browptr = new int[b_per_coldim+1];
		browptr[0] = 0;
		auto bcolinds = new int[bnnz];
		auto data = new T[bnnz*bm*bn];
		auto m = this->n;
		auto n = this->m;
		auto bm = this->bn;
		auto bn = this->bm;
		auto b_per_rowdim = this->b_per_coldim;
		auto b_per_coldim = this->b_per_rowdim;
		auto bnnz = this->bnnz;
		// copy-reorganize data in bmat.data
		// by scanning column by column of this to find each one of its blocks to rearrange in tbmat.data
		int t_block_offset = 0;
		for(int j=0;j<this->b_per_coldim;j++)
		{
			int i = 0; // number of blocks found so far in this j-th block column
			iter_block([this, &op, &j, &i, &t_block_offset, &data, &bcolinds, &browptr](int mat_row_id, int mat_col_id, int block_offset)
					{
					if(mat_col_id/this->bn == j)
					{
					DenseMatMap<T> block(data+t_block_offset*this->bn*this->bm, this->bn, this->bm);
					if(op == 'H')
						block = DenseMatMap<T>(this->data+block_offset*this->bm*this->bn, this->bm, this->bn).adjoint().eval();
					else if(op == 'T')
						block = DenseMatMap<T>(this->data+block_offset*this->bm*this->bn, this->bm, this->bn).adjoint().eval();
					bcolinds[t_block_offset++] = mat_row_id/this->bm;
					i++;
					}
					});
			browptr[j+1] = i;
		}
		for(int i=1;i<b_per_rowdim+1;i++)
			browptr[i] += browptr[i-1];
		free_bufs();
		this->m = m;
		this->n = n;
		this->bm = bm;
		this->bn = bn;
		this->b_per_rowdim = b_per_rowdim;
		this->b_per_coldim = b_per_coldim;
		this->bnnz = bnnz;
		this->browptr = browptr;
		this->bcolinds = bcolinds;
		this->data = data;
		return *this;
	}
	else
	{
		BSRMat<T, BlockStorageOrder> tbmat(*this);
		if(op == 'T')
			tbmat.transpose(/*inplace=*/true);
		else
			tbmat.adjoint(/*inplace=*/true);
		return tbmat;
	}
}

template<typename T, int BlockStorageOrder>
bool BSRMat<T, BlockStorageOrder>::contains_nan() const
{
	auto data_sz = bnnz*bm*bn;
	for(int i=0;i<data_sz;i++)
		if(std::isnan(std::real(data[i])))
			return true;
	return false;
}

template<typename T, int BlockStorageOrder>
const T& BSRMat<T, BlockStorageOrder>::operator()(unsigned int i, unsigned int j) const
{
	if(i < m && j < n)
	{
		// determine if (i,j) is inside a nz block
		unsigned int bi = i / bm; // the block row index of (i,j)
		// how many nz blocks in bi-th row 
		unsigned int bc = browptr[bi+1] - browptr[bi];
		if(!bc) return 0; // no nz blocks in the bi-th row
		// is j located in a nz block of bi-th row?
		unsigned int bj = j / bn; // block column index of (i,j)
		for(int k=browptr[bi];k<browptr[bi+1]; k++)
		{
			if(bcolinds[k] == bj)
				// (i,j) is in a nz block
				// return its value
				return data[k*bm*bn+j%bn*bm+i%bm];
		}
		// (i,j) is not in a nz block
		return zero;
	}
	throw std::runtime_error("BSRMat::operator() i or j is out of bounds.");
}

template<typename T, int BlockStorageOrder>
std::list<std::pair<int,int>> BSRMat<T, BlockStorageOrder>::nonzeros_indices() const
{
	std::list<std::pair<int,int>> nz_inds;
	iter_block([&nz_inds, this](int mat_row_id, int mat_col_id, int block_offset)
			{
				for(int i=0;i<bm;i++)
					for(int j=0;j<bn;j++)
					{
						if(data[block_offset*bm*bn+j*bm+i] != T(0)) // nz item
							nz_inds.push_back(std::make_pair(mat_row_id+i, mat_col_id+j));
					}
			});
	return nz_inds;
}

template<typename T, int BlockStorageOrder>
size_t BSRMat<T, BlockStorageOrder>::nbytes() const
{
	return bm*bn*bnnz*sizeof(T) /* data byte size */ + sizeof(int) * (bnnz + b_per_rowdim+1) /* indices byte size */;
}

template<typename T, int BlockStorageOrder>
size_t BSRMat<T, BlockStorageOrder>::nnz() const
{
	auto data_sz = bm*bn*bnnz;
	auto data_map = DenseMatMap<T>(data, data_sz, 1);
	auto nnz = 0;
	for(int i=0;i<data_sz;i++)
		if(data_map(i,0) != T(0))
			nnz++;
	return nnz;
}

template<typename T, int BlockStorageOrder>
Real<T> BSRMat<T, BlockStorageOrder>::density() const
{
	return nnz() / Real<T>(m*n);
}

template<typename T, int BlockStorageOrder>
Vec<T> BSRMat<T, BlockStorageOrder>::get_col(unsigned int col_id) const
{
	return Vec<T>(this->get_cols(col_id, 1));
}


template<typename T, int BlockStorageOrder>
SparseMat<T> BSRMat<T, BlockStorageOrder>::get_rows(unsigned int start_row_id, unsigned int num_rows) const
{
	if (start_row_id > m || start_row_id + num_rows > m)
		throw std::runtime_error("BSRMat::get_cols: matrix index overflow");
	SparseMat<T> smat(num_rows, n);
	typedef Eigen::Triplet<T,int> TRI;
	std::vector<TRI> triplets;
	iter_block([&triplets, this, &start_row_id, &num_rows](int mat_row_id, int mat_col_id, int block_offset)
			{
				for(int row_id=start_row_id;row_id<start_row_id+num_rows;row_id++)
				{
					if(mat_row_id <= row_id && row_id < mat_row_id+bm)
					{
							for(int j = 0; j < bn; j++)
							{
								int i = row_id%bm;
								triplets.push_back(Eigen::Triplet<T>(row_id-start_row_id, j+mat_col_id, (T)data[block_offset*bm*bn+j*bm+i]));
							}
					}
				}
			});
	smat.setFromTriplets(triplets.begin(), triplets.end());
	smat.makeCompressed();
	return smat;
}

template<typename T, int BlockStorageOrder>
SparseMat<T> BSRMat<T, BlockStorageOrder>::get_rows(const unsigned long int *row_ids, unsigned int num_rows) const
{
	for(int i=0;i<num_rows;i++)
		if (row_ids[i] > m)
			throw std::runtime_error("BSRMat::get_rows: matrix index overflow");
	SparseMat<T> smat(num_rows, n);
	typedef Eigen::Triplet<T,int> TRI;
	std::vector<TRI> triplets;
	iter_block([&triplets, this, &row_ids, &num_rows](int mat_row_id, int mat_col_id, int block_offset)
			{
			for(int i=0;i<num_rows;i++)
			{
				if(row_ids[i] >= mat_row_id && row_ids[i] < mat_row_id+bm)
				{
					for(int j=0; j < bn; j++)
					{
						int k = row_ids[i]%bm;
						triplets.push_back(Eigen::Triplet<T>(i, j+mat_col_id, (T)data[block_offset*bm*bn+j*bm+k]));
						}
				}
			}
			});
	smat.setFromTriplets(triplets.begin(), triplets.end());
	smat.makeCompressed();
	return smat;
}

template<typename T, int BlockStorageOrder>
SparseMat<T> BSRMat<T, BlockStorageOrder>::get_cols(unsigned int start_col_id, unsigned int num_cols) const
{
	if (start_col_id > n || start_col_id + num_cols > n)
		throw std::runtime_error("BSRMat::get_cols: matrix index overflow");
	SparseMat<T> smat(this->m, num_cols);
	typedef Eigen::Triplet<T,int> TRI;
	std::vector<TRI> triplets;
	iter_block([&triplets, this, &start_col_id, &num_cols](int mat_row_id, int mat_col_id, int block_offset)
			{
				for(unsigned int col_id=start_col_id;col_id < start_col_id+num_cols;col_id++)
					if(mat_col_id <= col_id &&  col_id < mat_col_id+bn)
					{
						unsigned int j = col_id%bn;
						for(int i=0; i < bm; i++)
							triplets.push_back(Eigen::Triplet<T>(i+mat_row_id, col_id-start_col_id, (T)data[block_offset*bm*bn+j*bm+i]));
					}
			});
	smat.setFromTriplets(triplets.begin(), triplets.end());
	smat.makeCompressed();
	return smat;
}

template<typename T, int BlockStorageOrder>
SparseMat<T> BSRMat<T, BlockStorageOrder>::get_cols(const unsigned long int *col_ids, unsigned int num_cols) const
{
	for(int j=0;j<num_cols;j++)
		if (col_ids[j] > n)
			throw std::runtime_error("bsrmat::get_cols: matrix index overflow");
	SparseMat<T> smat(this->m, num_cols);
	typedef Eigen::Triplet<T,int> TRI;
	std::vector<TRI> triplets;
	iter_block([&triplets, this, col_ids, &num_cols](int mat_row_id, int mat_col_id, int block_offset)
			{
			for(int j=0;j<num_cols;j++)
			{
				if(col_ids[j] >= mat_col_id && col_ids[j] < mat_col_id+bn)
				{
					unsigned int k = col_ids[j]%bn;
					for(int i=0; i < bm; i++)
					{
						triplets.push_back(Eigen::Triplet<T>(i+mat_row_id, j, (T) data[block_offset*bm*bn+k*bm+i]));
					}
				}
			}
			});
	smat.setFromTriplets(triplets.begin(), triplets.end());
	smat.makeCompressed();
	return smat;
}

template<typename T, int BlockStorageOrder>
BSRMat<T, BlockStorageOrder>::~BSRMat()
{
	free_bufs();
}


template<typename T, int BlockStorageOrder>
void BSRMat<T, BlockStorageOrder>::free_bufs()
{
	if(data != nullptr)
		delete[] data;
	if(browptr != nullptr)
		delete[] browptr;
	if(bcolinds != nullptr)
		delete[] bcolinds;
	browptr = bcolinds = nullptr;
	data = nullptr;
}
