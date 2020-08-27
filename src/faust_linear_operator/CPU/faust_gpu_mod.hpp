template <typename FPP>
int* Faust::FaustGPU<FPP>::gm_users = 0;

template <typename FPP>
void* Faust::FaustGPU<FPP>::gm_handle = nullptr;

template <typename FPP>
void* Faust::FaustGPU<FPP>::marr_funcs = nullptr;

template <typename FPP>
void* Faust::FaustGPU<FPP>::dsm_funcs = nullptr;

template <typename FPP>
void* Faust::FaustGPU<FPP>::gp_funcs = nullptr;


template <typename FPP>
std::map<void*,void*> Faust::FaustGPU<FPP>::cpu_gpu_map;

	template <typename FPP>
void Faust::FaustGPU<FPP>::check_gpu_mod_loaded()
{
	if(gm_handle == nullptr)
		throw std::runtime_error("Faust::enable_gpu_mod() must be called before any use of FaustGPU.");
}

	template <typename FPP>
void* Faust::FaustGPU<FPP>::init_gpu_mod(const std::string& libpath, const bool silent, void* gm_handle)
{
	if(Faust::FaustGPU<FPP>::gm_handle == nullptr)
		if(gm_handle == nullptr)
			Faust::FaustGPU<FPP>::gm_handle = gm_load_lib(libpath.c_str(), silent);
		else
			Faust::FaustGPU<FPP>::gm_handle = gm_handle;
	else
		if(! silent)
			std::cerr << "Warning: gm_lib is already loaded (can't reload it)." << endl;
	return Faust::FaustGPU<FPP>::gm_handle;
}

template<typename FPP>
bool Faust::FaustGPU<FPP>::is_cpu_mat_known(const Faust::MatGeneric<FPP,Cpu> *m)
{
	return cpu_gpu_map.find(const_cast<MatGeneric<FPP,Cpu>*>(m)) != cpu_gpu_map.end();
}

template<typename FPP>
bool Faust::FaustGPU<FPP>::are_cpu_mat_all_known(const std::vector<MatGeneric<FPP,Cpu>*> mats)
{
	for(auto m: mats)
		if(! is_cpu_mat_known(m)) return false;
	return true;
}

namespace Faust
{
	template<>
		void set_one<double>(double* scal)
		{
			*scal = 1.;
		}

	template<>
		void set_one<cuDoubleComplex>(cuDoubleComplex* scal)
		{
			scal->x = 1.;
			scal->y = 0;
		}
}

void* Faust::enable_gpu_mod(const char* libpath, const bool silent) //TODO: add backend argument: cuda (only available for now), opencl
{
	void* gm_handle = Faust::FaustGPU<double>::init_gpu_mod(std::string(libpath), silent, nullptr);
	if(gm_handle != nullptr)
		Faust::FaustGPU<std::complex<double>>::init_gpu_mod(std::string(libpath), silent, gm_handle);
	return gm_handle;
}

#include "faust_gpu_mod_double.hpp"
#include "faust_gpu_mod_complexdouble.hpp"
