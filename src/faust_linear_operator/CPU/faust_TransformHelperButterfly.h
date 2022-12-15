#ifndef __FAUST_TRANSFORM_HELPER_DFT__
#define __FAUST_TRANSFORM_HELPER_DFT__
#define IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
#include "faust_TransformHelper.h"

#include "faust_MatButterfly.h"
#include "faust_MatPerm.h"

#include <memory> // shared_ptr

namespace Faust
{
	template<typename FPP, FDevice DEV>
		class TransformHelperButterfly;

	template<typename FPP>
		class TransformHelperButterfly<FPP, Cpu> : public TransformHelper<FPP, Cpu>
		{
			using VecMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>;
			using DiagMat = Eigen::DiagonalMatrix<FPP, Eigen::Dynamic>;
			bool has_permutation;
			MatPerm<FPP, Cpu> perm;
			std::vector<std::shared_ptr<Faust::MatButterfly<FPP, Cpu>>> opt_factors;


			// private ctors
			TransformHelperButterfly<FPP, Cpu>(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
			TransformHelperButterfly(const TransformHelperButterfly<FPP,Cpu>* th, bool transpose, bool conjugate);


			public:
			std::string to_string() const;
			Vect<FPP, Cpu> multiply(const Vect<FPP, Cpu>& x);
			void multiply(const FPP* x, FPP* y);
			Vect<FPP,Cpu> multiply(const FPP* x);
			void multiply(const FPP* A, int A_ncols, FPP* C);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &A);
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A);
			static TransformHelper<FPP,Cpu>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP, Cpu>* optFaust(const TransformHelper<FPP, Cpu>* F);
			TransformHelper<FPP,Cpu>* transpose();
		};

}
#include "faust_TransformHelperButterfly.hpp"
#endif
