#ifndef __FAUST_GPU2CPU_H__
#define __FAUST_GPU2CPU_H__
#include "faust_cuda.h"


//modif AL AL
#include "faust_constant_gpu.h"
//#include "faust_MatDense_gpu.h"


// modif AL AL
template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class Vect;

#ifdef __COMPILE_SPMAT__
    template<typename FPP,Device DEVICE> class MatSparse;
    template<typename FPP,Device DEVICE> class Transform;
#endif

//! \brief faust_gpu2cpu copy an CPU vector/matrix into a GPU vector/matrix depend of the overloading function.
template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Vect<FPP,Cpu>& v, const Faust::Vect<FPP1,Gpu>& cu_v, cudaStream_t stream=0);
template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::MatDense<FPP,Cpu>& M, const Faust::MatDense<FPP1,Gpu>& cu_M, cudaStream_t stream=0);
#ifdef __COMPILE_SPMAT__
    template<typename FPP, typename FPP1>
    void faust_gpu2cpu(Faust::MatSparse<FPP,Cpu>& S, const Faust::MatSparse<FPP1,Gpu>& cu_S, cudaStream_t stream=0);
    template<typename FPP, typename FPP1>
    void faust_gpu2cpu(Faust::Transform<FPP,Cpu>& fcore, const Faust::Transform<FPP1,Gpu>& cu_fcore);
#endif


#include "faust_gpu2cpu.hpp"

#endif
