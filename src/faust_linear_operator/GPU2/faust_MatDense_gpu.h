#ifndef __FAUST_MATDENSE_GPU2__
#define __FAUST_MATDENSE_GPU2__
#ifdef USE_GPU_MOD
#include "faust_MatDense.h"
#include "faust_gpu_mod_utils.h"
namespace Faust
{

	template<typename FPP,FDevice DEVICE>
    class MatDense;

	template <typename FPP>
	void gemm(const MatDense<FPP, GPU2> &A, const MatDense<FPP, GPU2> &B, MatDense<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB);

	template<typename FPP>
		class MatDense<FPP, GPU2> : MatDense<FPP, Cpu>
		{
			friend void gemm<>(const MatDense<FPP, GPU2> &A, const MatDense<FPP, GPU2> &B, MatDense<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA, const char opB);

			public:
				MatDense(const faust_unsigned_int nbRow,
						const faust_unsigned_int nbCol,
						const FPP* data = nullptr,
						const bool no_alloc=false,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatDense(const MatDense<FPP,Cpu>& mat, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/);

				void operator=(const MatDense<FPP,GPU2> & A);
				void operator=(const MatDense<FPP,Cpu> & A);
				void operator=(const MatSparse<FPP,Cpu> & A);
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
				//  other = (*this) * other
				void multiply(MatDense<FPP, GPU2> &other, const char op_this='N');
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
				void operator*=(const MatDense<FPP, GPU2> &other);
				void operator*=(const MatDense<FPP, Cpu> &other);
//				void operator*=(MatSparse<FPP, Cpu> &other);
				void operator*=(const FPP& lambda);
				void scalarMultiply(const FPP& lambda);
				void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol);
				void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
				void setOnes();
				void setZeros();
				void setEyes();
				void transpose();
				void conjugate();
				void adjoint();
				void abs();
				Real<FPP> spectralNorm(const faust_unsigned_int nbr_iter_max, const float threshold);
				Real<FPP> norm();
				Real<FPP> normL1();
				void normalize();
				MatDense<FPP, GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr);
				void move(const int32_t dev_id=-1, const void* stream=nullptr);
				MatDense<FPP, Cpu> tocpu(const void* stream=nullptr);
				~MatDense<FPP, GPU2>();
			private:
				gm_DenseMat_t gpu_mat;
		};

};
#include "faust_MatDense_gpu_double.hpp"
#endif
#endif
