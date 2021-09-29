#ifndef __FAUST_MAT_DIAG__
#define __FAUST_MAT_DIAG__
#include "matio.h"
#include <Eigen/Core>
#include "faust_constant.h"
#include "faust_Vect.h"
#include "faust_MatGeneric.h"
#include <functional>
#include <list>
#include <utility>
#include  <algorithm>
#include <exception>
#include <complex>
#include "faust_linear_algebra.h"

namespace Faust
{
//	template<typename FPP, FDevice DEVICE>
//		class MatGeneric;

	template<typename FPP>
		class MatDiag : public MatGeneric<FPP,Cpu>
		{

			Eigen::DiagonalMatrix<FPP,Eigen::Dynamic> mat;
			bool isZeros; //TODO: init!

			public:
			MatDiag(faust_unsigned_int n): MatGeneric<FPP,Cpu>(n,n), isZeros(true) {}
			MatDiag(faust_unsigned_int nrows, faust_unsigned_int ncols): MatGeneric<FPP,Cpu>(nrows,ncols), isZeros(true) {}

			MatDiag(faust_unsigned_int n, const FPP* data): MatGeneric<FPP,Cpu>(n,n), isZeros(false)
			{
				Eigen::Matrix<FPP, Eigen::Dynamic, 1> v(n);
				memcpy(v.data(), data, sizeof(FPP)*n);
				mat = v.asDiagonal();
//				cout << mat.diagonal().size() << endl;
//				cout << mat.rows() << " " << mat.cols() << endl;
				isZeros = getNonZeros() == 0;
			}

#ifdef _MSC_VER
			MatDiag(faust_unsigned_int nrows, faust_unsigned_int ncols, const FPP* data): MatDiag<FPP>(min(nrows,ncols), data)
#else
			MatDiag(faust_unsigned_int nrows, faust_unsigned_int ncols, const FPP* data): MatDiag<FPP>(std::min(nrows,ncols), data)
#endif

			{
				this->dim1 = nrows;
				this->dim2 = ncols;
			}

			MatDiag(const MatDiag<FPP> & M) : MatDiag(M.dim1, M.dim2, M.getData())
			{
			}


			MatType getType() const { return Diag; }

			MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;

			void multiply(Vect<FPP,Cpu> & vec, char opThis='N') const;

			void multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			void multiply(MatSparse<FPP,Cpu> & M, char opThis) const { throw std::exception();}
			void multiplyRight(MatSparse<FPP,Cpu> const & M) { throw std::bad_function_call();}
			void transpose() { faust_unsigned_int tmp; tmp = this->dim1; this->dim1 = this->dim2; this->dim2 = tmp; }
			void conjugate(const bool eval = true) { mat = mat.diagonal().conjugate().asDiagonal(); }
			void adjoint() { conjugate(true); transpose(); }
			faust_unsigned_int getNonZeros() const { return mat.diagonal().nonZeros(); }

			size_t getNBytes() const;
			matvar_t* toMatIOVar(bool transpose, bool conjugate) const;
			Real<FPP> normL1(const bool transpose=false) const;
			typename Eigen::NumTraits<FPP>::Real norm() const;
			Real<FPP> normL1(faust_unsigned_int& col_id, const bool transpose) const;

			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;

			void operator*=(const FPP alpha);
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &v) const;
			void faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const;


			MatGeneric<FPP,Cpu>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatGeneric<FPP,Cpu>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const{if(i == j) return getData()[i]; else throw std::runtime_error("MatDiag::operator()(int,int): row and column indices must be equal.");}
			std::list<std::pair<int,int>> nonzeros_indices() const;
			//! \brief Returns all the features of the MatDense.
			std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;
			void Display() const;
			void setZeros();

			bool containsNaN();
			const FPP* getData() const { return mat.diagonal().data();};

		};
}

#include "faust_MatDiag.hpp"
#endif
