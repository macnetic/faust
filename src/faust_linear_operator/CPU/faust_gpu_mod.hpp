
template <typename FPP>
std::map<void*,void*> Faust::FaustGPU<FPP>::cpu_gpu_map;

template<typename FPP>
bool Faust::FaustGPU<FPP>::is_cpu_mat_known(const Faust::MatGeneric<FPP,Cpu> *m)
{
	return cpu_gpu_map.find(const_cast<MatGeneric<FPP,Cpu>*>(m)) != cpu_gpu_map.end();
}

template<typename FPP>
bool Faust::FaustGPU<FPP>::are_cpu_mat_all_known(const std::vector<MatGeneric<FPP,Cpu>*> mats)
{
	if(mats.size() == 0) return false;
	for(auto m: mats)
		if(! is_cpu_mat_known(m)) return false;
	return true;
}

namespace Faust
{
	//TODO: these functions have to move elsewhere more appropriate or simply disappear
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



#include "faust_gpu_mod_double.hpp"
#include "faust_gpu_mod_complexdouble.hpp"
