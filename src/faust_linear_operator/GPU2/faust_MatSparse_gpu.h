#ifndef __FAUST_MATSPARSE_GPU2__
#define __FAUST_MATSPARSE_GPU2__
#ifdef USE_GPU_MOD
#define __GM_LOADER__
#include "gm_interf.h"
#include "faust_constant.h"
namespace Faust
{

	template<typename FPP,FDevice DEVICE>
    class MatSparse;


	template<typename FPP>
		class MatSparse<FPP, GPU2>
		{

			public:
				/** \brief Inits from CPU buffers.
				 *
				 */
				MatSparse(const faust_unsigned_int nbRow,
						const faust_unsigned_int nbCol,
						const int32_t nnz = 0,
						const FPP* values = nullptr,
						const int32_t* rowptr = nullptr,
						const int32_t* colinds = nullptr,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatSparse<FPP, Cpu>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatDense<FPP, Cpu>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatSparse<FPP, GPU2>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				void operator=(const MatSparse<FPP, GPU2>& mat);

				void tocpu(MatSparse<FPP,Cpu> &sp_mat);
				Real<FPP> norm();
				void transpose();
				void adjoint();
				void conjugate();

				void resize(int32_t nnz, int32_t nrows, int32_t ncols);
				void setEyes();
				void setZeros();
				MatSparse<FPP, GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr);
				int32_t getNbRow();
				int32_t getNbCol();
				int32_t getNonZeros();
				~MatSparse();

			private:
				int32_t nbRow;
				int32_t nbCol;
				static void* dsm_funcs;
				static void* spm_funcs;
				static void* gp_funcs;
				gm_SparseMat_t gpu_mat;
		};


	template <typename FPP>
		void* Faust::MatSparse<FPP,GPU2>::dsm_funcs = nullptr;

	template <typename FPP>
		void* Faust::MatSparse<FPP,GPU2>::spm_funcs = nullptr;

	template <typename FPP>
		void* Faust::MatSparse<FPP,GPU2>::gp_funcs = nullptr;

};
#include "faust_MatSparse_gpu_double.hpp"
#endif
#endif
