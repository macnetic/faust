#ifndef __FAUST_MATGENERIC_GPU__
#define __FAUST_MATGENERIC_GPU__
#include "faust_constant.h"
#include "faust_Timer.h"
#include "faust_MatDense_gpu.h"
#include "faust_MatSparse_gpu.h"
#include "faust_MatPerm_gpu.h"
#include "faust_MatButterfly_gpu.h"

namespace Faust
{
	template<typename FPP, FDevice DEVICE> class Transform;
	template<typename FPP, FDevice DEVICE> class MatGeneric;
	template<typename FPP>
		class Transform<FPP,GPU2>;
	// The interest of this class is mostly to make Transform capable of storing generic matrix
	// TODO: this class should extends MatGeneric<FPP,Device>
	template<typename FPP>
		class MatGeneric<FPP, GPU2>
		{
			friend Transform<FPP,GPU2>; // needs access to get_gpu_mat_ptr
			virtual void set_gpu_mat_ptr(void*)=0;
			protected:
				bool is_identity;
				bool is_zeros;
			public:
				virtual void setZeros()=0;
				virtual size_t getNBytes() const=0;
				virtual MatType getType() const=0;
				virtual int32_t getNbRow() const=0;
				virtual int32_t getNbCol() const=0;
				virtual MatGeneric<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const=0;
				virtual MatGeneric<FPP,GPU2>* Clone(const bool isOptimize=false) const=0;
				virtual void* get_gpu_mat_ptr() const=0;
				virtual faust_unsigned_int getNonZeros() const=0;
				Real<FPP> density() const{return ((Real<FPP>) this->getNonZeros())/((float)this->getNbCol()*this->getNbRow());}
				virtual void transpose()=0;
				virtual void conjugate()=0;
				virtual void adjoint()=0;
				//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
				virtual MatGeneric<FPP,GPU2>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const=0;
				//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
				virtual MatGeneric<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const=0;
				//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
				virtual Faust::MatGeneric<FPP,GPU2>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const=0;
				//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
				virtual Faust::MatGeneric<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const=0;

				//! \brief Display the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
				virtual void Display() const;

				//! \brief Returns the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
				virtual std::string to_string(MatType type, const bool transpose=false, const bool displaying_small_mat_elts=false) const;
				//
				//! \brief Returns the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
				virtual std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;

				std::string to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity, MatType type) const;

			static std::string get_scalar_type_str();
				virtual Real<FPP> norm() const=0;

				virtual void multiply(MatDense<FPP,GPU2> &A, const char opThis) const =0;
				virtual void operator*=(const FPP& alpha) =0;

				// \brief Converts the matrix to a MatDense<FPP, DEVICE>.
				virtual MatDense<FPP, GPU2> to_dense() const=0;
				MatGeneric();

				virtual ~MatGeneric();
		};

//	template<typename FPP>
//	MatGeneric<FPP,GPU2>* optimize(MatDense<FPP,GPU2> const & M, MatSparse<FPP,GPU2> const & S);
//
}
#include "faust_MatGeneric_gpu.hpp"
#endif
