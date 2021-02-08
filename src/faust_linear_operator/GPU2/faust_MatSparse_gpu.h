#ifndef __FAUST_MATSPARSE_GPU2__
#define __FAUST_MATSPARSE_GPU2__
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatSparse.h"
namespace Faust
{
	template<typename FPP, FDevice DEVICE>
		class MatSparse;
	template<typename FPP, FDevice DEVICE>
		class MatDense;

	template<typename FPP>
		class MatSparse<FPP, GPU2> : public MatGeneric<FPP,GPU2>
		{

			friend Transform<FPP,GPU2>; // need to access to get_gpu_mat_ptr
			friend MatDense<FPP,GPU2>;
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
						const void* stream=nullptr,
						const bool nozero=false);

				MatSparse(const MatSparse<FPP, Cpu>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatDense<FPP, Cpu>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatSparse<FPP, GPU2>& mat,
						const int32_t dev_id=-1,
						const void* stream=nullptr);

				MatSparse(const MatDense<FPP,GPU2>& mat);

				MatSparse(MatSparse<FPP,GPU2> && mat);
				MatSparse<FPP,GPU2>& operator=(MatSparse<FPP,GPU2> && mat);

				MatSparse();

				MatSparse<FPP,GPU2>& operator=(const MatSparse<FPP, GPU2>& mat);
				void operator=(const MatSparse<FPP, Cpu>& mat);
				void operator=(const MatDense<FPP, GPU2>& mat);
				void operator*=(const FPP& alpha);
				void operator/=(const FPP& alpha);
				/**
				 * Subtracts a scalar value ONLY to nonzeros (zeros stay zeros).
				 */
				void operator-=(const FPP& alpha);
				/**
				 * Adds a scalar value ONLY to nonzeros (zeros stay zeros).
				 */
				void operator+=(const FPP& alpha);
				bool operator==(const MatSparse<FPP, GPU2>& mat) const;
				bool operator!=(const MatSparse<FPP, GPU2>& mat) const;

				void tocpu(MatSparse<FPP,Cpu> &sp_mat) const;
				void tocpu(int* row_ptr, int* col_ind, FPP* value_ptr, int* nrows=nullptr, int* ncols=nullptr, int* nnz=nullptr) const;
				Real<FPP> norm();
				void transpose();
				void adjoint();
				void conjugate();

				void resize(int32_t nnz, int32_t nrows, int32_t ncols);
				void setEyes();
				void setIdentity(int32_t dim);
				void setZeros();
				void set(int32_t nnz, int32_t nrows, int32_t ncols, FPP* values, int32_t* rowptr, int32_t* colids);
				void set(int32_t nnz, int32_t nrows, int32_t ncols, FPP* values, size_t* rowptr, size_t* colids);
				MatSparse<FPP, GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const;
				MatGeneric<FPP,GPU2>* Clone(const bool isOptimize=false) const;
				void move(const int32_t dev_id=-1, const void* stream=nullptr);
				int32_t getNbRow() const;
				int32_t getNbCol() const;
				size_t getNBytes() const;
				faust_unsigned_int getNonZeros() const;
				int32_t getDevice() const;
				void Display() const;
				std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;
				MatType getType() const;
				void multiply(Vect<FPP,GPU2>& vec, char opThis='N') const;
				~MatSparse();

			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatSparse<FPP,GPU2>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			MatSparse<FPP,GPU2>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;


			private:
				void* get_gpu_mat_ptr() const;
				void set_gpu_mat_ptr(void*);
				gm_SparseMat_t gpu_mat;
		};


};

#endif
#endif
