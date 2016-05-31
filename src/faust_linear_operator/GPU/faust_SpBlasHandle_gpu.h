#ifndef SPBLASHANDLE_GPU_H
#define SPBLASHANDLE_GPU_H
#include "faust_constant.h"
#include "cuda.h"
#include "cusparse.h"


namespace Faust
{

    template <Device DEVICE> class SpBlasHandle;

    template <>
    class SpBlasHandle<Gpu>
    {

    public :

        SpBlasHandle(){}
        SpBlasHandle(cusparseHandle_t cusparseHandle):m_cusparseHandle(cusparseHandle){}
        SpBlasHandle(const SpBlasHandle<Gpu> & spblashandleGPU): m_cusparseHandle(spblashandleGPU.m_cusparseHandle){}
        void operator=(SpBlasHandle<Gpu> const& spblashandleGPU){}//{m_cusparseHandle=spblashandleGPU.m_cusparseHandle;}
        void Init(){cusparseStatus_t cusparseStat=cusparseCreate(&m_cusparseHandle);}
        void Destroy(){cusparseDestroy(m_cusparseHandle);}
        const cusparseHandle_t & GetHandle() const {return m_cusparseHandle;}
        private :
        cusparseHandle_t m_cusparseHandle;


    };

}

#endif
