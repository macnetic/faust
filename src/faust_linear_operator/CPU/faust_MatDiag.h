#ifndef __FAUST_MAT_DIAG__
#define __FAUST_MAT_DIAG__
#include <Eigen/Core>
#include "faust_constant.h"
#include "faust_Vect.h"
#include "faust_MatGeneric.h"
using namespace Eigen;
using namespace std;


namespace Faust
{
//	template<typename FPP, Device DEVICE>
//		class MatGeneric;

	template<typename FPP>
		class MatDiag : public MatGeneric<FPP,Cpu>
		{

			DiagonalMatrix<FPP,Dynamic> mat;

			public:
			MatDiag(faust_unsigned_int n): MatGeneric<FPP,Cpu>(n,n) {}

			MatDiag(faust_unsigned_int n, FPP* data): MatGeneric<FPP,Cpu>(n,n)
			{
				Matrix<FPP, Dynamic, 1> v(n);
				memcpy(v.data(), data, sizeof(FPP)*n);
				mat = v.asDiagonal();
			}

			MatType getType() const { return Diag; }

			MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;

			void multiply(Vect<FPP,Cpu> & vec, char opThis='N') const;

			void multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			void transpose() { /* nothing to do */ }
			void conjugate() {  };//TODO:
			faust_unsigned_int getNonZeros() const { return 0; } //TODO

			matvar_t* toMatIOVar(bool transpose, bool conjugate) const;
			FPP normL1(const bool transpose) const;
			FPP norm() const;
			FPP normL1(faust_unsigned_int& col_id, const bool transpose) const;

			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;

			void operator*=(const FPP alpha);
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &v) const;
			void faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const;


			MatGeneric<FPP,Cpu>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatGeneric<FPP,Cpu>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;
			MatGeneric<FPP,Cpu>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
		};
#include "faust_MatDiag.hpp"
}
#endif
