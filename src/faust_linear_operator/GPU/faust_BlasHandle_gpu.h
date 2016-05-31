#ifndef BLASHANDLE_GPU_H
#define BLASHANDLE_GPU_H
#include "faust_constant.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"


namespace Faust
{

    template <Device DEVICE> class BlasHandle;


    template <>
    class BlasHandle<Gpu>
    {

    public :
        BlasHandle(cublasHandle_t cublasHandle):m_cublasHandle(cublasHandle){}
        BlasHandle(){}
        BlasHandle(const BlasHandle<Gpu> & blashandleGPU): m_cublasHandle(blashandleGPU.m_cublasHandle){}
        void operator=(BlasHandle<Gpu> const& blashandleGPU){m_cublasHandle=blashandleGPU.m_cublasHandle;}
        void Init(){cublasStatus_t cublasStat=cublasCreate(&m_cublasHandle);}
        void Destroy(){cublasDestroy(m_cublasHandle);}
        const cublasHandle_t & GetHandle() const{return m_cublasHandle;}

        private :
        cublasHandle_t m_cublasHandle;


    };

}

#endif
