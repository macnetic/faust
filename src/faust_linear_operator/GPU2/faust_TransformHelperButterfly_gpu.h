#ifndef __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#define __FAUST_TRANSFORM_HELPER_DFT_GPU2__
#include "faust_TransformHelper_gpu.h"

namespace Faust
{
	template<typename FPP, FDevice DEV>
		class TransformHelperButterfly;

	template<typename FPP>
	class ButterflyMat;

	template<typename FPP>
		class TransformHelperButterfly<FPP, GPU2> : public TransformHelper<FPP, GPU2>
		{
//			using VecMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, 1>>;
//			using DiagMat = Eigen::DiagonalMatrix<FPP, Eigen::Dynamic>;
//			FPP *perm_d_ptr;
//			DiagMat D;
//			std::vector<unsigned int> bitrev_perm;
//			std::vector<ButterflyMat<FPP>> opt_factors;
//
//
//			// private ctor
//			TransformHelperButterfly<FPP, GPU2>(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
//

			public:
			static TransformHelperButterfly<FPP,GPU2>* fourierFaust(unsigned int n, const bool norma=true) { throw std::runtime_error("Not yet implemented on GPU");};
//			Vect<FPP, GPU2> multiply(const Vect<FPP, GPU2>& x);
//			void multiply(const FPP* x, FPP* y);
//			Vect<FPP,GPU2> multiply(const FPP* x);
//			void multiply(const FPP* A, int A_ncols, FPP* C);
//			MatDense<FPP, GPU2> multiply(const MatDense<FPP,GPU2> &A);
//			MatDense<FPP, GPU2> multiply(const MatSparse<FPP,GPU2> &A);

		};

}
//#include "faust_TransformHelperButterfly_gpu.hpp" //TODO
#endif
