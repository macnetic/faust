#ifndef LINALGEBRA_GPU2_H
#define LINALGEBRA_GPU2_H

#include "faust_constant.h"
#include "faust_MatDense_gpu.h"
#include "faust_Vect_gpu.h"

namespace Faust
{

	// Computes alpha*typeA(A)*typeB(B)+ beta*C into C.
	template<typename FPP>
		void gemm(const MatDense<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP  alpha, const FPP  beta, char  opA, char  opB);

	// Computes alpha*opA(A)*b+ beta*c into c.
	template<typename FPP>
		void gemv(const MatDense<FPP, GPU2> &A, const Vect<FPP, GPU2> &b, Vect<FPP, GPU2> &c, const FPP& alpha, const FPP& beta, const char opA);

	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	template<typename FPP>
		void gemm_gen(const MatGeneric<FPP,GPU2> & A, const MatGeneric<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP  alpha, const FPP  beta, char  opA, char  opB);

	//	TODO: implements using MatSparse::multiply, warning: 'H' is not supported for opB (see gpu_mod / https://docs.nvidia.com/cuda/archive/9.2/cusparse/index.html cusparseTcsrmm2 for more details), so do a copy-conjugate manually beforehand
//	template<typename FPP>
//		void spgemm(const MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char opA, char opB);
//
//	template<typename FPP>
//		void spgemm(const MatDense<FPP,Cpu> & A,const MatSparse<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char opA, char opB);

}
#include "faust_linear_algebra_gpu.hpp"
#endif
