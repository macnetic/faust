#ifndef __FAUST_MATBSR_GPU2__
#define __FAUST_MATBSR_GPU2__
#ifdef USE_GPU_MOD
#include <complex>
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatBSR.h"

namespace Faust
{
	template<typename FPP, FDevice DEVICE>
		class MatSparse;
	template<typename FPP, FDevice DEVICE>
		class MatDense;
	template<typename FPP, FDevice DEVICE>
		class MatBSR;

	template<typename FPP>
		class MatBSR<FPP, GPU2> : public MatGeneric<FPP,GPU2>
		{


			gm_BSRMat_t gpu_mat;
			/** \brief Inits from CPU buffers.
			 *
			 */

			public:
			/*********** ctors **************/
			MatBSR(const faust_unsigned_int nrows,
					const faust_unsigned_int ncols,
					const faust_unsigned_int bnrows,
					const faust_unsigned_int bncols,
					const faust_unsigned_int bnnz,
					const FPP* bdata,
					const int32_t* browptr,
					const int32_t* bcolinds,
					const int32_t dev_id=-1,
					const void* stream=nullptr);

			MatBSR(const MatBSR<FPP, Cpu>& mat,
					const int32_t dev_id=-1,
					const void* stream=nullptr);

			MatBSR();

			/*********** MatGeneric member functions **************/
			void setZeros();
			size_t getNBytes() const;
			MatType getType() const;
			int32_t getNbRow() const;
			int32_t getNbCol() const;
			faust_unsigned_int getNonZeros() const;
			int32_t getDevice() const;
			MatGeneric<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;
			MatGeneric<FPP,GPU2>* Clone(const bool isOptimize=false) const;
			void* get_gpu_mat_ptr() const;
			void transpose();
			void conjugate();
			void adjoint();
			MatGeneric<FPP,GPU2>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatGeneric<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			Faust::MatGeneric<FPP,GPU2>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			Faust::MatGeneric<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;

			void Display() const;
			std::string to_string() const;
			Real<FPP> norm() const;

			virtual void set_gpu_mat_ptr(void*);
			/*********** own member functions **************/
			void multiply(MatDense<FPP,GPU2>& M, char opThis='N') const;
			void multiply(Vect<FPP,GPU2>& vec, char opThis='N') const;

			void tocpu(int32_t* browptr, int32_t* bcolinds, FPP* bdata, int32_t* nrows=nullptr, int32_t* ncols=nullptr, int32_t *bnrows=nullptr, int32_t *bncol=nullptr, int32_t* bnnz=nullptr) const;
			void tocpu(MatBSR<FPP, Cpu> &cpu_mat) const;
			MatSparse<FPP, GPU2> to_sparse() const;

		};
}
#include "faust_MatBSR_gpu.hpp"
#endif
#endif
