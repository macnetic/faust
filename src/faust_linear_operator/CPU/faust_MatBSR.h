#ifndef __BSR_MAT__
#define __BSR_MAT__
#include <Eigen/SparseCore>
#include "matio.h"
#include <Eigen/Dense>
#include <list>
#include <utility>
#include "faust_MatGeneric.h"

//template<typename FPP, FDevice DEVICE> class MatGeneric;

template<typename T>
using DenseMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using DenseMatMap = Eigen::Map<DenseMat<T>>;

template<typename T>
using SparseMat = Eigen::SparseMatrix<T,Eigen::RowMajor>;

#ifndef FAUST_CONSTANT_H
template<typename T>
using Real = typename Eigen::NumTraits<T>::Real;
#endif

template<typename T,int BlockStorageOrder=0> class BSRMat;

namespace Faust
{
	template<typename FPP>
		class MatBSR<FPP,Cpu> : public MatGeneric<FPP,Cpu>
		{
			BSRMat<FPP> bmat; // low-level BSRMat
			MatBSR() : MatGeneric<FPP, Cpu>() {}
			MatBSR(BSRMat<FPP>& mat);
			public:
			MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;
			void multiply(Vect<FPP,Cpu> & vec, char opThis) const;
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &v) const; // from LinearOperator
			void multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			void multiply(MatSparse<FPP, Cpu>& M, char opThis) const;
			void multiplyRight(MatSparse<FPP, Cpu> const& M);
			void faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const; // from LinearOperator
			void transpose();
			void conjugate(const bool eval=true);
			void adjoint();
			faust_unsigned_int getNonZeros()const;
			size_t getNBytes() const;
			MatType getType() const;
			void operator*=(const FPP alpha);
			std::string to_string(MatType type, const bool transpose=false, const bool displaying_small_mat_elts=false) const;
			std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;
			matvar_t* toMatIOVar(bool transpose, bool conjugate) const;
			Real<FPP> normL1(const bool transpose=false) const;
			Real<FPP> norm() const;
			Real<FPP> normL1(faust_unsigned_int& col_id, const bool transpose) const;
			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;
			MatGeneric<FPP,Cpu>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatGeneric<FPP,Cpu>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			std::list<std::pair<int,int>> nonzeros_indices() const;
			void setZeros();
			bool containsNaN();
			const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const ;
			MatDense<FPP, Cpu> to_dense() const;
			/**
			 * Returns the number of nonzeros blocks.
			 */
			size_t getNBlocks() const;
			size_t getNbBlocksPerDim(int dim_id) const;
			/**
			 * Returns number of rows of each nonzero block.
			 */
			size_t getNbBlockRow() const;
			/**
			 * Returns the number of columns of each nonzero block.
			 */
			size_t getNbBlockCol() const;
			virtual ~MatBSR();
			static MatBSR<FPP, Cpu>* randMat(int m, int n, int bm, int bn, int bnnz);
		};
}


//TODO: put all this class and member defs in Faust namespace after unit tests integration
template<typename T,int BlockStorageOrder>
class BSRMat
{

	friend class Faust::MatBSR<T, Cpu>;
	// \brief the data of nonzero blocks (size bnnz*bm*bn).
	T* data;
	// \brief column indices of the nonzero blocks (size bnnz).
	int* bcolinds;// TODO: should be unsigned
	// \brief the buffer of cumulative number of blocks per block row (size b_per_coldim+1).
	int* browptr;// TODO: should be unsigned
	// \brief the number of nonzero blocks
	int bnnz;// TODO: should be unsigned
	// \brief the number of matrix rows.
	int m; // TODO: should be unsigned
	// \brief the number of matrix columns.
	int n;// TODO: should be unsigned
	// \brief the number of block rows.
	int bm;// TODO: should be unsigned
	// \brief the number of block columns.
	int bn;// TODO: should be unsigned
	// \brief the number of blocks along the row dimension of the matrix (m/bm).
	int b_per_rowdim;// TODO: should be unsigned
	// \brief the number of blocks along the column dimension of the matrix (n/bn).
	int b_per_coldim;// TODO: should be unsigned

	// private default constructor
	BSRMat(): data(nullptr), bcolinds(nullptr), browptr(nullptr), bnnz(0), m(0), n(0), bm(0), bn(0), b_per_rowdim(0), b_per_coldim(0) {}
	public:
	/** Copy constructor */
	BSRMat(const BSRMat<T, BlockStorageOrder>& src_bmat);
	/** Move constructor */
	BSRMat(BSRMat<T, BlockStorageOrder>&& src_bmat);

	/** Copy operator */
	BSRMat<T, BlockStorageOrder>& operator=(const BSRMat<T, BlockStorageOrder>& src_bmat);
	/** Move operator */
	BSRMat<T, BlockStorageOrder>& operator=(BSRMat<T, BlockStorageOrder>&& src_bmat);


	/**
	 * Prints internal buffers of this: browptr, bcolinds, data.
	 */
	void print_bufs();
	/**
	 * Returns true if this contains NaN elements.
	 */
	bool contains_nan() const;
	/**
	 * Returns the (i,j) entry of this.
	 */
	T operator()(unsigned int i, unsigned int j) const;
	/**
	 * Returns the indices on the nonzeros of this.
	 */
	std::list<std::pair<int,int>> nonzeros_indices() const;
	/**
	 * Returns the number of bytes consumed by this.
	 */
	size_t nbytes() const;
	/**
	 * Returns the nnz of this.
	 */
	size_t nnz() const;
	/**
	 * Returns the density of this.
	 */
	Real<T> density() const;
	/**
	 * Returns one column of this as a vector.
	 *
	 * \param col_id: the column index.
	 */
	Vec<T> get_col(unsigned int col_id) const;
	/**
	 * Returns a sequence of num_rows of this into a SparseMat.
	 *
	 * \param start_row_id: the first row to return.
	 * \param num_rows: the number of rows returned.
	 */
	SparseMat<T> get_rows(unsigned int start_row_id, unsigned int num_rows) const;
	/**
	 * Returns a sequence of num_cols columns of this into a SparseMat.
	 *
	 * \param start_col_id: the first column to return.
	 * \param num_cols: the number of columns returned.
	 */
	SparseMat<T> get_cols(unsigned int start_col_id, unsigned int num_cols) const;
	/**
	 * Returns a sequence of rows of this into a SparseMat.
	 *
	 * \param row_ids: the indices of the rows in the desired order.
	 * \param num_rows: the size of row_ids.
	 */
	SparseMat<T> get_rows(const unsigned long int* row_ids, unsigned int num_rows) const;
	/**
	 * Returns a sequence of columns of this into a SparseMat.
	 *
	 * \param col_ids: the indices of the columns in the desired order.
	 * \param num_cols: the size of col_ids.
	 */
	SparseMat<T> get_cols(const unsigned long int* col_ids, unsigned int num_cols) const;
	/**
	 * Converts this to a DenseMat.
	 */
	DenseMat<T> to_dense() const;
	/**
	 * Converts this to a SparseMat.
	 */
	SparseMat<T> to_sparse() const;
	/**
	 * Multiplies this by a DenseMat and returns the result as a DenseMat.
	 */
	DenseMat<T> mul(const DenseMat<T>&) const;
	/**
	 * Multiplies this by a SparseMat and returns the result as a DenseMat.
	 */
	DenseMat<T> mul(const SparseMat<T>&) const;
	/**
	 * Multiplies this by a vector (passed directly as a buffer) and returns the result as a Vec.
	 */
	Vec<T> mul(const T* vec_data) const;
	/**
	 * Multiplies this by a BSRMat and returns the result as DenseMat.
	 */
	DenseMat<T> mul(const BSRMat<T, BlockStorageOrder>&) const;
	/**
	 * Multiplies this by a scalar (in place).
	 */
	void mul(const T& scal) const;
	/**
	 * Iterates on the matrix blocks, row by row.
	 * For each block the argument function is called with the row and column indices of the upper left corner of the block and the block_offset in data and bcolinds.
	 */
	void iter_block(std::function<void(int /* mat_row_id */, int /* mat_col_id */, int /* block_offset */)> op) const;
	/**
	 * Frobenius norm.
	 */
	Real<T> norm() const;

	/**
	 * Inf-norm.
	 */
	Real<T> normInf() const;

	/**
	 * 1-norm
	 */
	Real<T> normL1() const;
	/**
	 * Transpose
	 */
	BSRMat<T, BlockStorageOrder> transpose(const bool inplace=false);
	BSRMat<T, BlockStorageOrder> conjugate(const bool inplace=false);
	/**
	 * \param m: matrix number of rows.
	 * \param n: matrix number of columns.
	 * \param bm: number of rows of blocks (must divide m evenly).
	 * \param bn: number of columns of blocks (must divide n evenly).
	 * \param bnnz: number of blocks that contain nonzeros among the m*n/bm/bn blocks.
	 * \param min_bdensity: the inner minimal density of each block.
	 * \param max_bdensity: the inner maximal density of each block.
	 */
	static BSRMat<T, BlockStorageOrder> rand(int m, int n, int bm, int bn, int bnnz);

	void free_bufs();
	~BSRMat();
};
#include "faust_MatBSR.hpp"
#endif
