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

	// Computes alpha*opA(A)*b+ beta*c into C.
	template<typename FPP>
		void gemv(const MatDense<FPP, GPU2> &A, const Vect<FPP, GPU2> &b, Vect<FPP, GPU2> &C, const FPP& alpha, const FPP& beta, const char opA);

	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	template<typename FPP>
		void gemm_gen(const MatGeneric<FPP,GPU2> & A, const MatGeneric<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP  alpha, const FPP  beta, char  opA, char  opB);

	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	template<typename FPP>
		void spgemm(const MatSparse<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB);

	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	// \param impl_meth: in any case this function rely on previous spgemm prototype, if impl_meth is 1 then transpose/transconjugate is used to avoid converting A and B to another type of matrix, otherwise (impl_meth is any other value) A is converted to a MatSparse and B to a MatDense
	template<typename FPP>
		void spgemm(const MatDense<FPP,GPU2> & A,const MatSparse<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB, int impl_meth = 1);
	//
	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	template<typename FPP>
		void bsrgemm(const MatBSR<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB);

	// Computes alpha*opA(A)*opB(B)+ beta*C into C.
	// \param impl_meth: in any case this function rely on previous bsrgemm prototype, if impl_meth is 1 then transpose/transconjugate is used to avoid converting A and B to another type of matrix, otherwise (impl_meth is any other value) A is converted to a MatSparse and B to a MatDense
	template<typename FPP>
		void bsrgemm(const MatDense<FPP,GPU2> & A,const MatBSR<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB, int impl_meth = 1);
}
#include "faust_linear_algebra_gpu.hpp"
#endif
