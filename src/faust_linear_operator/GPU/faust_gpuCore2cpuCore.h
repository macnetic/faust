#ifndef __FAUST_GPUCORE_2_CPUCORE_H__
#define __FAUST_GPUCORE_2_CPUCORE_H__



template<typename FPP,Device DEVICE> class Transform;

template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Transform<FPP,Cpu>& fcore, const Faust::Transform<FPP,Gpu>& cu_fcore);

#include "faust_gpuCore2cpuCore.hpp"
#endif
