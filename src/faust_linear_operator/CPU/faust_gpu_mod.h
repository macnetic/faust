#ifndef __FAUST_GPU__
#define __FAUST_GPU__
#define __GM_LOADER__
#include "gm_interf.h"
#include <exception>
#include <iostream>
#include <vector>
#include "faust_MatGeneric.h"

namespace Faust
{
	template<typename FPP>
	class FaustGPU
	{

		static void *gm_handle;
		static int *gm_users;
		static void* marr_funcs; //void because we don't know FPP yet and templates aren't available through shared lib interface (extern C, no name mangling)
		static void* dsm_funcs;
		gm_MatArray_t gpu_mat_arr;
		size_t size;

		public:
		FaustGPU(std::vector<Faust::MatGeneric<FPP,Cpu>*>&);
		~FaustGPU();

		MatDense<FPP, Cpu> get_product();
		MatDense<FPP, Cpu> multiply(const Faust::MatGeneric<FPP,Cpu>*);

		static void* init_gpu_mod(const std::string& libpath = "libgm.so", const bool silent = false, void* gm_handle = nullptr);
		static void check_gpu_mod_loaded();
	};

	template<typename T>
		void set_one(T* scal);

	void* enable_gpu_mod(const char* libpath= "libgm.so", const bool silent = false);

}

#include "faust_gpu_mod.hpp"

#endif
