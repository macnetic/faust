#ifndef __MAT_BUTTERFLY_GPU__
#define __MAT_BUTTERFLY_GPU__
#include "faust_MatSparse.h"

namespace Faust
{

	template<typename FPP, FDevice DEVICE>
		class MatButterfly;

	template<typename FPP>
		class MatButterfly<FPP, GPU2> : public MatGeneric<FPP,GPU2>
		{


			Vect<FPP, GPU2> d1;
			Vect<FPP, GPU2> d2;
			int* subdiag_ids;
			int level;

			Vect<FPP, GPU2> d2t;
			bool is_transp;

			public:
			//TODO: move defs in hpp
			MatButterfly() : level(-1), subdiag_ids(nullptr), d1(), d2(), d2t(), is_transp(false)
			{
			}

			MatButterfly(const MatButterfly<FPP, GPU2>& bmat)
			{
				*this = bmat;
			}

			MatButterfly& operator=(const MatButterfly<FPP, GPU2>& bmat)
			{
				this->d1 = bmat.d1;
				this->d2 = bmat.d2;
				this->d2t = bmat.d2t;
				this->level = bmat.level;
				this->subdiag_ids = new int[d1.size()];
				std::copy(bmat.subdiag_ids, bmat.subdiag_ids+d1.size(), this->subdiag_ids);
				if(bmat.is_transp)
					transpose();
				return *this;
			}

			MatButterfly(const MatSparse<FPP, GPU2> &factor, int level)
			{
				MatSparse<FPP, Cpu> Scpu;
				factor.tocpu(Scpu);
				MatButterfly<FPP, GPU2> this_(Scpu, level);
				*this = this_;
			//TODO: do it without passing through CPU mem. and move the def in hpp
			}


			MatButterfly(const MatButterfly<FPP, Cpu>& bmat) : MatButterfly(bmat.toMatSparse(), bmat.getLevel())
			{
				//TODO/ without conversion to MatSparse
			}

			MatDense<FPP, GPU2> multiply(const FPP* x);
			MatDense<FPP, GPU2> multiply(const FPP* A, int A_ncols);
			MatDense<FPP, GPU2> multiply(MatDense<FPP,GPU2> &A);
			void multiply(MatDense<FPP,GPU2> &A, MatDense<FPP, Cpu> & out);
			const Vect<FPP, GPU2>& getD1() {return d1;};
			const Vect<FPP, GPU2>& getD2() {return d2;};
			const int getLevel() const {return level;}

			MatButterfly(const MatSparse<FPP, Cpu> &factor, int level);
			void setZeros();
			size_t getNBytes() const;
			MatType getType() const {return Butterfly;} //TODO: move def in hpp
			int32_t getNbRow() const {return d1.size();} //TODO: move def in hpp
			int32_t getNbCol() const {return d1.size();} //TODO: move def in hpp
			faust_unsigned_int getNonZeros() const;
			MatButterfly<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;
			MatButterfly<FPP,GPU2>* Clone(const bool isOptimize=false) const;
			void transpose();
			void init_transpose();
			void conjugate();
			void adjoint();
			void* get_gpu_mat_ptr() const { throw std::runtime_error("get_gpu_mat_ptr doesn't make sense for a MatButterfly");}
			void set_gpu_mat_ptr(void*) { throw std::runtime_error("set_gpu_mat_ptr doesn't make sense for a MatButterfly");}
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;

			void Display() const;
			Real<FPP> norm() const;
			void multiply(MatDense<FPP, GPU2> &other, const char op_this);
			void multiply(MatSparse<FPP, GPU2> &other, const char op_this);
			MatSparse<FPP, GPU2> toMatSparse() const;
			~MatButterfly();
		};

}
#include "faust_MatButterfly_gpu.hpp"
#endif
