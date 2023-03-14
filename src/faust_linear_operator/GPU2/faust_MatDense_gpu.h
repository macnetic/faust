#ifndef __FAUST_MATDENSE_GPU2__
#define __FAUST_MATDENSE_GPU2__
#ifdef USE_GPU_MOD
#include <complex>
#include <cstdint>
#include "faust_MatDense.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_Vect_gpu.h"
#include "faust_MatSparse_gpu.h"
#include "faust_MatBSR_gpu.h"
#include "faust_gpu_mod_utils.h"
namespace Faust
{


	template<typename FPP, FDevice DEVICE>
		class MatDense;
	//TODO: move in a specific module
	template<typename FPP>
		void butterfly_diag_prod(MatDense<FPP, GPU2>& X, const Vect<FPP, GPU2>& d1, const Vect<FPP, GPU2>& d2, const int* ids);

	/***
	 * \brief Computes the SVD of a batch of matrices As.
	 *
	 * \param As: a batch of matrices horizontally concatenated.
	 * \param batch_sz: the number of matrices to compute the SVD of.
	 * \param Us: receives the left singulars vectors of each matrix of As (the singulars vectors of each matrix are concatenated horizontally).
	 * \param Vs: same as Us for the right singular vectors.
	 * \param Ss: revices the singular values for each matrix of As (the columns in Ss corresponds to the matrices of As in the same order).
	 */
	template<typename FPP>
		void batched_svd(MatDense<FPP, GPU2>& As, const uint32_t batch_sz, MatDense<FPP, GPU2>& Us, MatDense<FPP, GPU2>& Vs, MatDense<Real<FPP>, GPU2>& Ss, const uint32_t rank = 0);


	/***
	 * \brief Similar to the previous prototype except that the results is transferred into CPU matrices and only a number of rank singular vectors/values is copied.
	 *
	 */
	template<typename FPP>
		void batched_svd(MatDense<FPP, Cpu>& As, const uint32_t batch_sz, MatDense<FPP, Cpu>& Us, MatDense<FPP, Cpu>& Vs, MatDense<Real<FPP>, Cpu>& Ss, const uint32_t rank = 0);


	template<typename FPP>
		class MatDense<FPP, GPU2> : public MatGeneric<FPP,GPU2>
		{
			friend Transform<FPP,GPU2>; // need to access to get_gpu_mat_ptr
			friend MatSparse<FPP,GPU2>;
			friend MatBSR<FPP,GPU2>;
			friend MatDense<std::complex<double>,GPU2>; // TODO limit to real function
			friend void butterfly_diag_prod<>(MatDense<FPP, GPU2>& X, const Vect<FPP, GPU2>& d1, const Vect<FPP, GPU2>& d2, const int* ids);
			friend void batched_svd<>(MatDense<FPP, GPU2>& As, const uint32_t batch_sz, MatDense<FPP, GPU2>& Us, MatDense<FPP, GPU2>& Vs, MatDense<Real<FPP>, GPU2>& Ss, const uint32_t rank /*= 0*/);
			friend void batched_svd<>(MatDense<std::complex<Real<FPP>>, GPU2>& As, const uint32_t batch_sz, MatDense<std::complex<Real<FPP>>, GPU2>& Us, MatDense<std::complex<Real<FPP>>, GPU2>& Vs, MatDense<Real<FPP>, GPU2>& Ss, const uint32_t rank /*= 0*/);
			friend void batched_svd<>(MatDense<FPP, Cpu>& As, const uint32_t batch_sz, MatDense<FPP, Cpu>& Us, MatDense<FPP, Cpu>& Vs, MatDense<Real<FPP>, Cpu>& Ss, const uint32_t rank /*= 0*/);
			friend void batched_svd<>(MatDense<std::complex<Real<FPP>>, Cpu>& As, const uint32_t batch_sz, MatDense<std::complex<Real<FPP>>, Cpu>& Us, MatDense<std::complex<Real<FPP>>, Cpu>& Vs, MatDense<Real<std::complex<Real<FPP>>>, Cpu>& Ss, const uint32_t rank /*= 0*/);
//			friend void gemm<>(const MatDense<FPP, GPU2> &A, const MatDense<FPP, GPU2> &B, MatDense<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB);
//
//			friend void gemv<>(const MatDense<FPP, GPU2> &A, const Vect<FPP, GPU2> &B, Vect<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB);

			public:
				MatDense(const faust_unsigned_int nbRow,
						const faust_unsigned_int nbCol,
						const FPP* cpu_data = nullptr,
						const bool no_alloc=false,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatDense();
				MatDense(const MatDense<FPP,Cpu>& mat, const int32_t dev_id=-1, const void* stream=nullptr);

				MatDense(MatDense<FPP,GPU2> && mat);
				MatDense(const MatDense<FPP,GPU2> & mat);
				MatDense(const MatSparse<FPP,GPU2> & mat);
				~MatDense();

				MatDense<FPP,GPU2>& operator=(MatDense<FPP,GPU2> && mat);
				MatDense<FPP,GPU2>& operator=(const MatDense<FPP,GPU2> & A);
				void operator=(const MatDense<FPP,Cpu> & A);
				void operator=(const MatSparse<FPP,Cpu> & A);
				void operator=(const MatSparse<FPP,GPU2> & A);
				// *this = *this + A
				void add(const MatDense<FPP,Cpu> & A);
				void add(const MatDense<FPP,GPU2> & A);
				void add(const MatSparse<FPP,Cpu> & A);
				void operator+=(const MatDense<FPP,GPU2> & A);
				void operator+=(const MatDense<FPP,Cpu> & A);
				void operator+=(const MatSparse<FPP,Cpu> & A);
				//TODO: add(const MatSparse<FPP,GPU2>)
				void sub(const MatDense<FPP,Cpu> & A);
				void sub(const MatDense<FPP,GPU2> & A);
				void sub(const MatSparse<FPP,Cpu> & A);
				void operator-=(const MatDense<FPP,GPU2> & A);
				void operator-=(const MatDense<FPP,Cpu> & A);
				void operator-=(const MatSparse<FPP,Cpu> & A);
				//TODO: sub(const MatSparse<FPP,GPU2>)
				//TODO: void add(MatSparse<FPP,GPU2> const& A);
				// vec = this * vec
				Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu> &vec);
				/**
				 * *this = *this * vec (element-wise multiplication)
				 * possible broadcasting
				 * \param ids: if not nullptr then *this rows are multiplied in this order to vec.
				 */
				void eltwise_mul(const Vect<FPP, GPU2> &vec, const int *ids=nullptr);
				//  TODO: other shouldn't be const if it is the output
				//  other = (*this) * other
				void multiply(MatDense<FPP, GPU2> &other, const char op_this='N') const;
				//  other = (*this) * other
				void multiply(MatDense<FPP, Cpu> &other, const char op_this='N');
//				void multiply(MatSparse<FPP, Cpu> &other, MatDense<FPP, GPU2>& output, const char op_this='N');
				void multiply(const MatSparse<FPP, Cpu> &other, MatDense<FPP, Cpu>& output, const char op_this='N');
				void multiply(const MatSparse<FPP, Cpu> &other, MatDense<FPP, GPU2>& output, const char op_this='N');
				//! \brief Replace (this) by (this) * A
				void multiplyRight(const MatDense<FPP, Cpu>& A);
				void multiplyRight(const MatDense<FPP, GPU2>& A);
				//! \brief Replace (this) by S * (this)
				void multiplyLeft(const MatSparse<FPP, Cpu>& S, const char transS='N');
				void multiplyLeft(const MatSparse<FPP, GPU2>& S, const char transS='N');
				void multiply(const Vect<FPP, GPU2>& vec, Vect<FPP, GPU2>& out_vec) const;
				//! \brief compute MatDense-vector multiplication
				//! \param vec : the vector
				//! \param opThis : character	
				//! vec = (*this) * vec if opThis='N'
				// vec = (*this)' * vec if opThis='T'
				void multiply(Vect<FPP,GPU2> & vec,const char opThis) const
				{gemv((*this),vec,vec,(FPP) 1.0, (FPP) 0.0,opThis);}

				void operator*=(const MatDense<FPP, GPU2> &other);
				void operator*=(const MatDense<FPP, Cpu> &other);
//				void operator*=(MatSparse<FPP, Cpu> &other);
				void operator*=(const FPP& lambda);
				void scalarMultiply(const FPP& lambda);
				void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol);
				void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
				void setOnes();
				void setZeros();
				void setRand();
				void setEyes();
				void setData(const FPP* data, int32_t nrows, int32_t ncols);
				void transpose();
				void conjugate();
				void adjoint();
				void abs();
				Real<FPP> spectralNorm(const faust_unsigned_int nbr_iter_max, const float threshold);
				Real<FPP> norm() const;
				Real<FPP> normL1() const;
				void normalize();
				int32_t getDevice() const;
				FPP trace() const;
				MatDense<FPP, GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;

				MatGeneric<FPP,GPU2>* Clone(const bool isOptimize=false) const;
				void move(const int32_t dev_id=-1, const void* stream=nullptr);
				void tocpu(FPP* cpu_buffer, const void* stream/*=nullptr*/) const;
				MatDense<FPP, Cpu> tocpu(const void* stream=nullptr) const;
				void tocpu(MatDense<FPP, Cpu> &cpu_mat, const void* stream=nullptr) const;
				MatType getType() const;
				size_t getNBytes() const;
				int32_t getNbRow() const;
				int32_t getNbCol() const;
				faust_unsigned_int getNonZeros() const;
				MatDense<FPP,GPU2>* get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows) const;
				MatDense<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int n) const;
				MatDense<FPP,GPU2>* get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const;
				MatDense<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int n) const;
				void copyBuf(FPP* dst_cpu_buf, const void* stream=nullptr) const;
				/** \brief copy a block (defined by offset and size) of this matrix into a CPU buffer
				 *
				 * \param offset: defines the offset of the block to copy (in number of scalar elements).
				 * \param size: number of scalar elements to copy from this matrix to the CPU buffer, starting from the offset.
				 */
				void copyBufBlock(FPP* dst_cpu_buf, uint32_t offset = 0, int32_t size = -1, const void* stream=nullptr) const;
				bool isReal() const;
				void prox_sp(int32_t k, bool normalized=false, bool pos=false) const;
				void prox_spcol(int32_t k, bool normalized=false, bool pos=false) const;
				void prox_splin(int32_t k, bool normalized=false, bool pos=false) const;
		  void prox_triu_sp(int32_t k, bool normalized=false, bool pos=false) const;
		  void prox_tril_sp(int32_t k, bool normalized=false, bool pos=false) const;
				void real(MatDense<Real<FPP>, GPU2>& real_mat) const;
				template<typename FPP2>
					MatDense<Real<FPP2>, GPU2> to_real() const;

				MatDense<FPP, GPU2> to_dense() const;
				static void gemm(const MatDense<FPP, GPU2> &A, const MatDense<FPP, GPU2> &B, MatDense<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB);
				static void gemv(const MatDense<FPP, GPU2> &A, const Vect<FPP, GPU2> &B, Vect<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB='N');

			protected:
				gm_DenseMat_t gpu_mat;
				void* get_gpu_mat_ptr() const;
				void set_gpu_mat_ptr(void*);
		};


}

#include "faust_MatDense_gpu.hpp"

#endif
#endif
