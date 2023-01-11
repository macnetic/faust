#ifndef __MAT_BUTTERFLY_GPU__
#define __MAT_BUTTERFLY_GPU__
#include "faust_MatSparse.h"

namespace Faust
{

	template<typename FPP, FDevice DEVICE>
		class MatButterfly;

	template<typename FPP>
		class MatButterfly<FPP, GPU2>/* : public MatGeneric<FPP,GPU2>*/
		{


			Vect<FPP, GPU2> d1;
			Vect<FPP, GPU2> d2;
			int* subdiag_ids;
			int level;

			Vect<FPP, GPU2> d2t;
			bool is_transp;

			public:
			MatButterfly(const MatSparse<FPP, Cpu> &factor, int level);
			void setZeros();
			size_t getNBytes() const;
			MatType getType() const {return Butterfly;} //TODO: move def in hpp
			int32_t getNbRow() const {return d1.size();} //TODO: move def in hpp
			int32_t getNbCol() const {return d1.size();} //TODO: move def in hpp
			MatButterfly<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;
			MatButterfly<FPP,GPU2>* Clone(const bool isOptimize=false) const;
			/*  void* get_gpu_mat_ptr() const;
			  faust_unsigned_int getNonZeros() const;
			  void transpose();
			  void conjugate();
			  void adjoint();*/
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;

			/*void Display() const;*/
			Real<FPP> norm() const;
			void multiply(MatDense<FPP, GPU2> &other, const char op_this);
			void multiply(MatSparse<FPP, GPU2> &other, const char op_this);
			MatSparse<FPP, GPU2> toMatSparse() const;
		};

}
#include "faust_MatButterfly_gpu.hpp"
#endif
