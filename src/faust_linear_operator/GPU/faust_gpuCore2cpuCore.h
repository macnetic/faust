#ifndef __FAUST_GPUCORE_2_CPUCORE_H__
#define __FAUST_GPUCORE_2_CPUCORE_H__

#include "faust_Transform.h"
#include "faust_Transform_gpu.h"
#include "faust_gpu2cpu.h"

template<typename FPP,Device DEVICE> class Transform;

template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Transform<FPP,Cpu>& fcore, const Faust::Transform<FPP,Gpu>& cu_fcore);

#include "faust_gpuCore2cpuCore.hpp"
#endif
