#ifndef __FAUST_CU_VEC_HPP__
#define __FAUST_CU_VEC_HPP__


#include "faust_Vect.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "faust_exception.h"
#include "faust_reduce_gpu.h"
#include "faust_cuda.h"
#include "kernels.h"
#include "faust_gpu2cpu.h"


using namespace std;

template <typename FPP>
bool Faust::Vect<FPP,Gpu>::equality(Faust::Vect<FPP,Gpu> const &x, FPP precision) const
{
 	Faust::Vect<FPP,Cpu> x1;
    Faust::Vect<FPP,Cpu> this1;
	faust_gpu2cpu(x1,x);
	faust_gpu2cpu(this1,(*this));
    return x1.equality(this1,precision);

}

template <typename FPP>
Faust::Vect<FPP,Gpu>::Vect(const FPP *data_,const faust_unsigned_int dim_, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    if(dataFromGPU)
        copyFromDevice(data_, dim_, dstDevice, srcDevice, stream);
    else
        copyFromHost(data_, dim_, dstDevice, stream);
}

template <typename FPP>
Faust::Vect<FPP,Gpu>::Vect(const Faust::Vect<FPP,Gpu>& cu_v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    copyFromDevice(cu_v.data, cu_v.dim, dstDevice, cu_v.device, stream);
}

template <typename FPP>
Faust::Vect<FPP,Gpu>::Vect(const Faust::Vect<FPP,Cpu>& v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    copyFromHost(v.getData(), v.size(), dstDevice, stream);
}

template <typename FPP>
Faust::Vect<FPP,Gpu>::Vect(const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    resize(dim_, dstDevice);
}


template <typename FPP>
const char * Faust::Vect<FPP,Gpu>::m_className = "Faust::Vect<FPP,Gpu>::";

template <typename FPP>
void Faust::Vect<FPP,Gpu>::_create(const faust_unsigned_int dim_, int device_)
{
    if(dim_<0)
        handleError(m_className, "_create : incorrect dimensions");


    if(dim_==0)
        data = NULL;

    else
    {
        int currentGPU;
        faust_cudaGetDevice(&currentGPU);
        faust_cudaSetDevice(device_);

        if (data == NULL)
            faust_cudaMalloc((void**)&data, (dim_)*sizeof(FPP));
        else
            handleError(m_className, "_create : data has already been allocated on GPU");


        faust_cudaSetDevice(currentGPU);
    }
    dim = dim_;
    device = device_;

}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::_clear()
{
    if(dim>0)
    {
        int currentGPU;
        faust_cudaGetDevice(&currentGPU);
        faust_cudaSetDevice(device);

        if (data != NULL)
            faust_cudaFree(data);
        else
            handleError(m_className, "_clear : data has already been deleted on GPU");

        faust_cudaSetDevice(currentGPU);
    }
    else if(data!=NULL)
        handleError(m_className, "_clear : data of vector shoould be NULL");
    dim = 0;
    data = NULL;
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::resize(const faust_unsigned_int dim_, int device_)
{
    if(dim_!=dim || (device_!=device))
    {
        _clear();
        _create(dim_, device_);
    }
    else
    {
        dim = dim_;
        device = device_;
    }
}


template <typename FPP>
void Faust::Vect<FPP,Gpu>::copyFromHost(const FPP *data_, const faust_unsigned_int dim_,  int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    if(data_==NULL || dim_<=0)
        handleError(m_className,"copyFromHost : NULL data pointer or incorrect value of dimension");
    resize(dim_, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(data, data_, dim_*sizeof(FPP), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::copyFromDevice(const FPP *data_, const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    if(data_==NULL || dim_<=0)
        handleError(m_className,"copyFromDevice : NULL data pointer or incorrect value of dimension");
    resize(dim_, dstDevice);
    faust_cudaMemcpyPeerAsync(data, dstDevice, data_, srcDevice, dim_*sizeof(FPP), stream);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::copyToHost(FPP *data_, const faust_unsigned_int dim_, cudaStream_t stream/*=0*/)const
{
    if(data==NULL || dim_<=0)
        handleError(m_className,"copyToHost : NULL data pointer or incorrect value of dimension");
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(data_, data, dim_*sizeof(FPP), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::copyToDevice(FPP *data_, const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
    if(data==NULL || dim_<=0)
        handleError(m_className,"copyToDevice : NULL data pointer or incorrect value of dimension");
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);
    faust_cudaMemcpyPeerAsync(data_, dstDevice, data, device, dim_*sizeof(FPP), stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::init(const Faust::Vect<FPP,Gpu>& cu_v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{

    if(cu_v.dim == 0 || cu_v.data==NULL )
    {
        // display warning as matrix is empty
        resize(cu_v.dim);
        data = NULL ;
        device = dstDevice;
        return;
    }
    copyFromDevice(cu_v.data, cu_v.dim, dstDevice, cu_v.device, stream);
}


template <typename FPP>
void Faust::Vect<FPP,Gpu>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    if(device == dstDevice)
        return;
    if(dim==0 || data==NULL)
    {
        device = dstDevice;
        return;
    }

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    FPP* data_tmp;
    faust_cudaMalloc((void**)&data_tmp, dim*sizeof(FPP));
    faust_cudaMemcpyPeerAsync(data_tmp, dstDevice, data, device, dim*sizeof(FPP), stream);
    faust_cudaFree(data);
    data = data_tmp;

    faust_cudaSetDevice(currentGPU);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator=(const Faust::Vect<FPP,Gpu>& cu_v)
{
    init(cu_v, cu_v.device);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator=(const Faust::Vect<FPP,Cpu>& v)
{
    if(v.size() == 0)
    {
        resize(0);
        data = NULL;
    }
    else
        copyFromHost(v.getData(), v.size(), device);
}

template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::max()const
{
    if(dim>0)
    {
        int currentGPU;
        faust_cudaGetDevice(&currentGPU);
        faust_cudaSetDevice(device);
        FPP val = faust_cu_max(data, dim);
        faust_cudaSetDevice(currentGPU);
        return val;
    }
    else
    {
        // display warning as vector is empty
        return (FPP)0.0;
    }
}

template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::min()const
{
    if(dim>0)
    {
        int currentGPU;
        faust_cudaGetDevice(&currentGPU);
        faust_cudaSetDevice(device);
        FPP val = faust_cu_min(data, dim);
        faust_cudaSetDevice(currentGPU);
        return val;
    }
    else
    {
        // display warning as vector is empty
        return (FPP)0.0;
    }
}

template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::sum()const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val = faust_cu_sum(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (FPP)0.0;
   }
}

template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::norm() const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val = faust_cu_norm(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (FPP)0.0;
   }
}

template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::dot(const Faust::Vect<FPP,Gpu>& cu_v, Faust::BlasHandle<Gpu> blasHandle) const
{
//return dot(*this, cu_v, blasHandle);
   if(size() != cu_v.size())
      handleError("LinAlgebra_cu","dot : the two vectors don't have the same size");
   if(size() > 0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      const Faust::Vect<FPP,Gpu>* cu_v_ptr = &cu_v;
      if(cu_v.getDevice() != device)
         cu_v_ptr = new Faust::Vect<FPP,Gpu>(cu_v, device);

      FPP result=0.0;

      faust_cu_dot(blasHandle.GetHandle(), size(), getData(), 1, cu_v_ptr->getData(), 1, &result);

      if(cu_v.getDevice() != device)
         delete cu_v_ptr;
      faust_cudaSetDevice(currentGPU);
      return result;
   }
   else
   {
      // display warning as vector is empty
      return (FPP)0.0;
   }


}

template <typename FPP>
bool Faust::Vect<FPP,Gpu>::operator==(const Faust::Vect<FPP,Gpu>& cu_v)const
{
   if(data!=NULL && data==cu_v.data && device== cu_v.device)
   {
      if(dim!=cu_v.dim)
         handleError(m_className,"operator== : same data and device but different dimensions");
      return true;
   }
   else
      return false;
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator+=(const Faust::Vect<FPP,Gpu>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(m_className,"operator+= : vector dimensions must agree");
   }
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   Faust::Vect<FPP,Gpu> vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_add(data, vec_tmp.data, dim);
   }
   else
      kernel_add(data, cu_v.data, dim);

   faust_cudaSetDevice(currentGPU);
}
template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator-=(const Faust::Vect<FPP,Gpu>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(m_className,"operator-= : vector dimensions must agree");
   }
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   Faust::Vect<FPP,Gpu> vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_sub(data, vec_tmp.data, dim);
   }
   else
      kernel_sub(data, cu_v.data, dim);
   faust_cudaSetDevice(currentGPU);
}
template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator*=(const Faust::Vect<FPP,Gpu>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(m_className,"operator*= : vector dimensions must agree");
   }
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);
   Faust::Vect<FPP,Gpu> vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_mult(data, vec_tmp.data, dim);
   }
   else
      kernel_mult(data, cu_v.data, dim);
   faust_cudaSetDevice(currentGPU);
}
template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator/=(const Faust::Vect<FPP,Gpu>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(m_className,"operator/= : vector dimensions must agree");
   }
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   Faust::Vect<FPP,Gpu> vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_div(data, vec_tmp.data, dim);
   }
   else
      kernel_div(data, cu_v.data, dim);

   faust_cudaSetDevice(currentGPU);
}


template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator+=(const FPP alpha)
{
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_add_const( data, alpha, dim);

      faust_cudaSetDevice(currentGPU);
   }
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator-=(const FPP alpha)
{
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_sub_const( data, alpha, dim);

      faust_cudaSetDevice(currentGPU);
   }
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator*=(const FPP alpha)
{
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_mult_const( data, alpha, dim);

      faust_cudaSetDevice(currentGPU);
   }
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::operator/=(const FPP alpha)
{
   if(dim == 0)
   {
      // display warning as vector is empty
      return;
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      kernel_div_const( data, alpha, dim);

      faust_cudaSetDevice(currentGPU);
   }
}



template <typename FPP>
void Faust::Vect<FPP,Gpu>::setValues(const FPP val)
{
   if (dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      kernel_memset(data, val, dim);
      faust_cudaSetDevice(currentGPU);
   }
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::setValues(const FPP val, const faust_unsigned_int length)
{
   resize(length);
   setValues(val);
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::setZeros()
{setValues(0.0);}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::setZeros(const faust_unsigned_int length)
{setValues(0.0, length);}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::setOnes()
{setValues(1.0);}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::setOnes(const faust_unsigned_int length)
{setValues(1.0, length);}


template <typename FPP>
FPP Faust::Vect<FPP,Gpu>::mean_relative_error(const Faust::Vect<FPP,Gpu>& v_ref)const
{
   if(v_ref.size() != size())
     handleError(m_className,"relative_error : sizes are different");

   if(dim == 0)
   {
      // display warning as vector is empty
      return (FPP)0.0;
   }
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   Faust::Vect<FPP,Gpu> tmp(size());
   kernel_relative_error(tmp.getData(), v_ref.getData(), data, size());

   faust_cudaSetDevice(currentGPU);
   return tmp.mean();
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::Display()const
{
   Faust::Vect<FPP,Cpu> v;
   faust_gpu2cpu(v,*this);
   v.Display();
}

template <typename FPP>
void Faust::Vect<FPP,Gpu>::print_file(const char* filename)const
{
   Faust::Vect<FPP,Cpu> v;
   faust_gpu2cpu(v,*this);
   v.print_file(filename);
}

#endif
