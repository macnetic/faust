#ifndef FAUST_MATPERM_H
#define FAUST_MATPERM_H

#ifdef USE_PYTHONIC
#include "numpy/_numpyconfig.h"
#include "ButFactor_matmul.hpp"
#include <pythonic/include/numpy/array.hpp>
#include <pythonic/numpy/array.hpp>
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/include/types/ndarray.hpp"

using namespace pythonic;

// Helper to create a float 1D array from a pointer
	template <typename T>
types::ndarray<T, types::pshape<long>> arrayFromBuf1D(T* fPtr, long size)
{
	auto shape = types::pshape<long>(size);
	return types::ndarray<T, types::pshape<long>>(fPtr,shape,types::ownership::external);
}
#endif


#include "faust_constant.h"
#include <Eigen/Core>
#include <memory> // shared_ptr

#define diag_conj(D) DiagMat(D.diagonal().conjugate()) // D is a DiagMat

namespace Faust
{
	template<typename FPP, FDevice DEVICE> class MatPerm;

	template<typename FPP>
		class MatPerm<FPP,Cpu> : public MatGeneric<FPP,Cpu>
		{

			using VecMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>;
			using DiagMat = Eigen::DiagonalMatrix<FPP, Eigen::Dynamic>;
			DiagMat D;
			DiagMat DT;
			std::vector<int> perm_ids;
			std::vector<int> perm_ids_T;
#ifdef USE_PYTHONIC
			long *perm_ids_ptr;
#endif
			bool is_transp;
			FPP zero;


			public:
			MatPerm<FPP,Cpu>(const MatSparse<FPP, Cpu> &factor);

			MatPerm(); // default ctor
			MatPerm(const MatPerm& src); // copy ctor
			MatPerm<FPP, Cpu>& operator=(const MatPerm& src); // assignment operator

			void init_transpose();
			void Display() const;

			faust_unsigned_int getNbRow() const { return D.rows();}; // a permutation matrix is square
			faust_unsigned_int getNbCol() const { return D.cols();};

			void multiply(const FPP* x, FPP* y, size_t size, bool transpose = false, bool conjugate  = false);
			void multiply(const FPP* A, int A_ncols, FPP* C, size_t size, bool transpose = false, bool conjugate  = false);

			const DiagMat& getD() {return D;}; //TODO/ move to .hpp
			const std::vector<int>& get_perm_ids() {return perm_ids;}

			MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;
			void multiply(Vect<FPP,Cpu> & vec, char opThis='N') const;
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &v) const; // from LinearOperator
			void multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			void multiply(MatSparse<FPP, Cpu>& M, char opThis) const;
			void multiplyRight(MatSparse<FPP, Cpu> const& M) ;
			void transpose();
			void conjugate(const bool eval=true);
			void adjoint();
			faust_unsigned_int getNonZeros()const;
			size_t getNBytes() const;
			MatType getType() const;
			void operator*=(const FPP alpha);
			matvar_t* toMatIOVar(bool transpose, bool conjugate, const char *var_name=nullptr) const;
			Real<FPP> normL1(const bool transpose=false) const;
			Real<FPP> norm() const;
			Real<FPP> normL1(faust_unsigned_int& col_id, const bool transpose) const;
			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;
			MatGeneric<FPP,Cpu>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatGeneric<FPP,Cpu>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			std::list<std::pair<int,int>> nonzeros_indices(const double& tol=0) const;
			void setZeros();
			bool containsNaN()const;
			const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const;

			MatSparse<FPP, Cpu> toMatSparse() const;

			void faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const; //from LinearOperator

			/**
			 * \brief Returns true if the MatSparse S is permutation (in which case it is eligible to a MatPerm conversion).
			 *
			 * If the matrix S is well-formed the cost is O(nnz^2) in the worst case.
			 */
			static bool isPerm(const MatSparse<FPP, Cpu> &S, bool verify_ones=true);

			MatDense<FPP, Cpu> to_dense() const; //TODO: rename toMatDense
		};

}

#include "faust_MatPerm.hpp"
#endif
