#ifndef __MAT_PERM_GPU__
#define __MAT_PERM_GPU__
#include "faust_MatSparse.h"

namespace Faust
{

	template<typename FPP, FDevice DEVICE>
		class MatPerm;

	template<typename FPP>
		class MatPerm<FPP, GPU2> : public MatGeneric<FPP,GPU2>
		{


			Vect<FPP, GPU2> d;
			int* perm_ids;

			Vect<FPP, GPU2> dt;
			int* perm_ids_T;
			bool is_transp;

			public:
			//TODO: move defs in hpp
			MatPerm() : perm_ids(nullptr), perm_ids_T(nullptr), d(), dt(), is_transp(false)
			{
			}

			MatPerm(const MatPerm<FPP, GPU2>& bmat) : MatPerm()
			{
				*this = bmat;
			}

			MatPerm& operator=(const MatPerm<FPP, GPU2>& bmat)
			{
				this->d = bmat.d;
				this->dt = bmat.dt;
				this->perm_ids = new int[d.size()];
				std::copy(bmat.perm_ids, bmat.perm_ids+d.size(), this->perm_ids);
				if(bmat.is_transp)
					transpose();
				return *this;
			}


			MatPerm(const MatSparse<FPP, GPU2> &factor) : MatPerm()
			{
				MatSparse<FPP, Cpu> Scpu;
				factor.tocpu(Scpu);
				MatPerm<FPP, GPU2> this_(Scpu);
				*this = this_;
			//TODO: do it without passing through CPU mem. and move the def in hpp
			}


			MatPerm(const MatPerm<FPP, Cpu>& bmat) : MatPerm(bmat.toMatSparse())
			{
				//TODO/ without conversion to MatSparse
			}

			MatDense<FPP, GPU2> multiply(const FPP* x);
			MatDense<FPP, GPU2> multiply(const FPP* A, int A_ncols);
			MatDense<FPP, GPU2> multiply(MatDense<FPP,GPU2> &A);
			void multiply(MatDense<FPP,GPU2> &A, MatDense<FPP, Cpu> & out);
			const Vect<FPP, GPU2>& getD() {return d;};

			MatPerm(const MatSparse<FPP, Cpu> &factor);
			void setZeros();
			size_t getNBytes() const;
			MatType getType() const {return Perm;} //TODO: move def in hpp
			int32_t getNbRow() const {return d.size();} //TODO: move def in hpp
			int32_t getNbCol() const {return d.size();} //TODO: move def in hpp
			faust_unsigned_int getNonZeros() const;
			MatPerm<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;
			MatPerm<FPP,GPU2>* Clone(const bool isOptimize=false) const;
			void transpose();
			void init_transpose();
			void conjugate();
			void adjoint();
			void* get_gpu_mat_ptr() const { throw std::runtime_error("get_gpu_mat_ptr doesn't make sense for a MatPerm");}
			void set_gpu_mat_ptr(void*) { throw std::runtime_error("set_gpu_mat_ptr doesn't make sense for a MatPerm");}
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
			void multiply(MatDense<FPP, GPU2> &other, const char op_this) const;
			void multiply(MatSparse<FPP, GPU2> &other, const char op_this);
			MatSparse<FPP, GPU2> toMatSparse() const;

			static bool isPerm(const MatSparse<FPP, GPU2> &S, bool verify_ones=true)
			{
				//TODO: do it without copy in CPU mem and move def in hpp
				MatSparse<FPP, Cpu> Scpu;
				S.tocpu(Scpu);
				return MatPerm<FPP, Cpu>::isPerm(Scpu, verify_ones);
			}

			~MatPerm();
		};

}
#include "faust_MatPerm_gpu.hpp"
#endif
