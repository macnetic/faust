#ifndef __FAUST_MATGENERIC_GPU__
#define __FAUST_MATGENERIC_GPU__
#include "faust_constant.h"
#include "faust_Timer.h"
#include "faust_MatDense_gpu.h"
#include "faust_MatSparse_gpu.h"

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

			virtual void Display() const=0;
			virtual Real<FPP> norm() const=0;
			MatGeneric();

			virtual ~MatGeneric();
		};

//	template<typename FPP>
//	MatGeneric<FPP,GPU2>* optimize(MatDense<FPP,GPU2> const & M, MatSparse<FPP,GPU2> const & S);
//
}
#include "faust_MatGeneric_gpu.hpp"
#endif
