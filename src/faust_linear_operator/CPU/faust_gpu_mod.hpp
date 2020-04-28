template <typename FPP>
int* Faust::FaustGPU<FPP>::gm_users = 0;

template <typename FPP>
void* Faust::FaustGPU<FPP>::gm_handle = nullptr;

template <typename FPP>
void* Faust::FaustGPU<FPP>::marr_funcs = nullptr;

template <typename FPP>
void* Faust::FaustGPU<FPP>::dsm_funcs = nullptr;



	template <typename FPP>
void Faust::FaustGPU<FPP>::check_gpu_mod_loaded()
{
	if(gm_handle == nullptr)
		throw std::runtime_error("FaustGPU::init_gpu_mod() must be called before any use of FaustGPU.");
}

	template <typename FPP>
void* Faust::FaustGPU<FPP>::init_gpu_mod(const std::string& libpath, void* gm_handle)
{
	if(Faust::FaustGPU<FPP>::gm_handle == nullptr)
		if(gm_handle == nullptr)
			Faust::FaustGPU<FPP>::gm_handle = gm_load_lib(libpath.c_str());
		else
			Faust::FaustGPU<FPP>::gm_handle = gm_handle;
	else
		std::cerr << "Warning: gm_lib is already loaded (can't reload it)." << endl;
	return Faust::FaustGPU<FPP>::gm_handle;
}

namespace Faust
{
	template<>
		void set_one<double>(double* scal)
		{
			*scal = 1.d;
		}

	template<>
		void set_one<cuDoubleComplex>(cuDoubleComplex* scal)
		{
			scal->x = 1.d;
			scal->y = 0;
		}
}

void* Faust::enable_gpu_mod(const char* libpath) //TODO: add backend argument: cuda (only available for now), opencl
{
	void* gm_handle = Faust::FaustGPU<double>::init_gpu_mod(std::string(libpath));
	if(gm_handle != nullptr)
		Faust::FaustGPU<std::complex<double>>::init_gpu_mod(std::string(libpath), gm_handle);
	return gm_handle;
}

#include "faust_gpu_mod_double.hpp"
#include "faust_gpu_mod_complexdouble.hpp"
