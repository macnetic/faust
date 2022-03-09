#ifndef __FAUST_GPU_MOD_UTILS__
#define __FAUST_GPU_MOD_UTILS__
#include <cstdint>
#include <cstdlib>
#include "gm_interf_types.h"
#ifndef __GM_LOADER__
// some headers contain function definitions (that's a hack to avoid writing again the code for loading the lib functions in user code)
// only the utils module must contain these definitions to load the functions in memory (otherwise it'll end up with a multiple definitions compiling error)
#include "gm_cuComplex.h"
#include "gm_interf_float.h"
#include "gm_interf_double.h"
#include "gm_interf_cuComplex.h"
#include "gm_interf_cuDoubleComplex.h"
#include "gm_interf_GenPurposeFuncs.h"
#endif
#include <iostream>
#include <string>
#include <complex>

using namespace std;

namespace Faust
{
	class GPUModHandler
	{
		static GPUModHandler* singleton; // no need to create/initialize more than one instance of the handler for a program

		// gm library descriptor
		void* gm_handle;

		// General functions (with no specific scalar type)
		gm_GenPurposeFunc* gp_funcs_;

		// MatArray funcs for all types
		gm_MatArrayFunc_double* marr_funcs_double;
		gm_MatArrayFunc_float* marr_funcs_float;
		gm_MatArrayFunc_cuComplex* marr_funcs_cuComplex;
		gm_MatArrayFunc_cuDoubleComplex* marr_funcs_cuDoubleComplex;

		// DenseMat funcs for all scalar types
		gm_DenseMatFunc_double* dsm_funcs_double;
		gm_DenseMatFunc_float* dsm_funcs_float;
		gm_DenseMatFunc_cuComplex* dsm_funcs_cuComplex;
		gm_DenseMatFunc_cuDoubleComplex* dsm_funcs_cuDoubleComplex;


		// SparseMat funcs for all scalar types
		gm_SparseMatFunc_double* spm_funcs_double;
		gm_SparseMatFunc_float* spm_funcs_float;
		gm_SparseMatFunc_cuComplex* spm_funcs_cuComplex;
		gm_SparseMatFunc_cuDoubleComplex* spm_funcs_cuDoubleComplex;

		// BSRMat funcs for all scalar types
		gm_BSRMatFunc_double* bsr_funcs_double;
		gm_BSRMatFunc_float* bsr_funcs_float;
		gm_BSRMatFunc_cuComplex* bsr_funcs_cuComplex;
		gm_BSRMatFunc_cuDoubleComplex* bsr_funcs_cuDoubleComplex;

		GPUModHandler();
		void load_gm_functions();
		void* init_gpu_mod(const string& libpath, bool silent, void* gm_handle);

		public:
		~GPUModHandler();
		static GPUModHandler* get_singleton(const bool silent=true);
		void* enable_gpu_mod(const string& libpath, bool silent); //TODO: add backend argument: cuda (only available for now), opencl
		void* enable_gpu_mod(const char* libpath, bool silent); //TODO: delete when python and matlab wrappers will use string
		void check_gpu_mod_loaded() const;
		bool is_gpu_mod_loaded() const;

		// functions to access the good scalar type gpu_mod functions specifiying the type by argument value (templates are not possible here because of the library C interface )
		// functions for sparse matrix (typically used in MatSparse<FPP,GPU2>)
		gm_SparseMatFunc_double* spm_funcs(const double &d) const;
		gm_SparseMatFunc_float* spm_funcs(const float &d) const;
		gm_SparseMatFunc_cuComplex* spm_funcs(const complex<float> &c) const;
		gm_SparseMatFunc_cuDoubleComplex* spm_funcs(const complex<double> &c) const;

		// functions for dense matrix (typically used in MatDense<FPP,GPU2>)
		gm_DenseMatFunc_double* dsm_funcs(const double &d) const;
		gm_DenseMatFunc_float* dsm_funcs(const float &d) const;
		gm_DenseMatFunc_cuComplex* dsm_funcs(const complex<float> &c) const;
		gm_DenseMatFunc_cuDoubleComplex* dsm_funcs(const complex<double> &c) const;

		// functions for sequence of matrices (aka a Faust)
		gm_MatArrayFunc_double* marr_funcs(const double &d) const;
		gm_MatArrayFunc_float* marr_funcs(const float &d) const;
		gm_MatArrayFunc_cuComplex* marr_funcs(const complex<float> &c) const;
		gm_MatArrayFunc_cuDoubleComplex* marr_funcs(const complex<double> &c) const;

		// functions for sparse matrix (typically used in MatBSR<FPP,GPU2>)
		gm_BSRMatFunc_double* bsr_funcs(const double &d) const;
		gm_BSRMatFunc_float* bsr_funcs(const float &d) const;
		gm_BSRMatFunc_cuComplex* bsr_funcs(const complex<float> &c) const;
		gm_BSRMatFunc_cuDoubleComplex* bsr_funcs(const complex<double> &c) const;

		gm_GenPurposeFunc* gp_funcs() const;
	};

	// easy access wrapper for the Faust namespace
	void* enable_gpu_mod(const char* libpath = "libgm.so", const bool silent=true);
	bool is_gpu_mod_enabled();
	void char2gm_Op(const char& c, gm_Op & op);
}
#define load_all_mat_funcs(type) \
			marr_funcs_##type = new gm_MatArrayFunc_##type();\
			dsm_funcs_##type = new gm_DenseMatFunc_##type();\
			spm_funcs_##type = new gm_SparseMatFunc_##type();\
			bsr_funcs_##type = new gm_BSRMatFunc_##type();\
			load_marr_funcs_##type(gm_handle, marr_funcs_##type);\
			load_dsm_funcs_##type(gm_handle, dsm_funcs_##type);\
			load_spm_funcs_##type(gm_handle, spm_funcs_##type);\
			load_bsr_funcs_##type(gm_handle, bsr_funcs_##type)

#define delete_mat_funcs(type) \
			delete marr_funcs_##type; \
			delete dsm_funcs_##type; \
			delete spm_funcs_##type; \
			delete bsr_funcs_##type
#endif
