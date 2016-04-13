#ifndef __FAUST_CU_VEC_HPP__
#define __FAUST_CU_VEC_HPP__


//#include "faust_cu_vec.h"
#include "faust_vec.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "faust_exception.h"
#include "faust_cu_reduce.h"
#include "faust_cuda.h"
#include "kernels.h"
#include "faust_cu2faust.h"

using namespace std;

template <typename faust_real>
bool faust_cu_vec<faust_real>::equality(faust_cu_vec<faust_real> const &x, faust_real precision) const
{
 	faust_vec<faust_real> x1;
	faust_vec<faust_real> this1;
	faust_cu2faust(x1,x);
	faust_cu2faust(this1,(*this));
	return x1.equality(this1,precision);	
	

}

template <typename faust_real>
 faust_cu_vec<faust_real>::faust_cu_vec(const faust_real *data_,const faust_unsigned_int dim_, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    if(dataFromGPU)
        copyFromDevice(data_, dim_, dstDevice, srcDevice, stream); 
    else
        copyFromHost(data_, dim_, dstDevice, stream); 
}

template <typename faust_real>
faust_cu_vec<faust_real>::faust_cu_vec(const faust_cu_vec<faust_real>& cu_v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
        copyFromDevice(cu_v.data, cu_v.dim, dstDevice, cu_v.device, stream);
}

template <typename faust_real>
faust_cu_vec<faust_real>::faust_cu_vec(const faust_vec<faust_real>& v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
        copyFromHost(v.getData(), v.size(), dstDevice, stream); 
}

template <typename faust_real>
faust_cu_vec<faust_real>::faust_cu_vec(const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim(0), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    resize(dim_, dstDevice);
}



template <typename faust_real>
const char * faust_cu_vec<faust_real>::class_name = "faust_cu_vec<faust_real>::";

template <typename faust_real>
void faust_cu_vec<faust_real>::_create(const faust_unsigned_int dim_, int device_)
{
   if(dim_<0)
       handleError(class_name, "_create : incorrect dimensions");


   if(dim_==0)
      data = NULL;

   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device_);

      if (data == NULL)
         faust_cudaMalloc((void**)&data, (dim_)*sizeof(faust_real));
      else
         handleError(class_name, "_create : data has already been allocated on GPU");
         
      
      faust_cudaSetDevice(currentGPU);
   }
   dim = dim_;
   device = device_;

}

template <typename faust_real>
void faust_cu_vec<faust_real>::_clear()
{
   if(dim>0)    
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      
      if (data != NULL)
         faust_cudaFree(data);
      else
         handleError(class_name, "_clear : data has already been deleted on GPU");
         
      faust_cudaSetDevice(currentGPU);
   }  
   else if(data!=NULL)
      handleError(class_name, "_clear : data of vector shoould be NULL");
   dim = 0;
   data = NULL;
}

template <typename faust_real>
void faust_cu_vec<faust_real>::resize(const faust_unsigned_int dim_, int device_)
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


template <typename faust_real>
void faust_cu_vec<faust_real>::copyFromHost(const faust_real *data_, const faust_unsigned_int dim_,  int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    if(data_==NULL || dim_<=0)
       handleError(class_name,"copyFromHost : NULL data pointer or incorrect value of dimension");
    resize(dim_, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(data, data_, dim_*sizeof(faust_real), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::copyFromDevice(const faust_real *data_, const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    if(data_==NULL || dim_<=0)
       handleError(class_name,"copyFromDevice : NULL data pointer or incorrect value of dimension");
    resize(dim_, dstDevice);
    faust_cudaMemcpyPeerAsync(data, dstDevice, data_, srcDevice, dim_*sizeof(faust_real), stream);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::copyToHost(faust_real *data_, const faust_unsigned_int dim_, cudaStream_t stream/*=0*/)const
{
    if(data==NULL || dim_<=0)
       handleError(class_name,"copyToHost : NULL data pointer or incorrect value of dimension");
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(data_, data, dim_*sizeof(faust_real), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::copyToDevice(faust_real *data_, const faust_unsigned_int dim_, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
    if(data==NULL || dim_<=0)
       handleError(class_name,"copyToDevice : NULL data pointer or incorrect value of dimension");
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);
    faust_cudaMemcpyPeerAsync(data_, dstDevice, data, device, dim_*sizeof(faust_real), stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::init(const faust_cu_vec<faust_real>& cu_v, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
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


template <typename faust_real>
void faust_cu_vec<faust_real>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
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

   faust_real* data_tmp;
   faust_cudaMalloc((void**)&data_tmp, dim*sizeof(faust_real));
   faust_cudaMemcpyPeerAsync(data_tmp, dstDevice, data, device, dim*sizeof(faust_real), stream);
   faust_cudaFree(data);
   data = data_tmp;

   faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::operator=(const faust_cu_vec<faust_real>& cu_v)
{

     init(cu_v, cu_v.device);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::operator=(const faust_vec<faust_real>& v)
{
      if(v.size() == 0)
      {
         resize(0);
         data = NULL;         
      }
      else
         copyFromHost(v.getData(), v.size(), device);
}

template <typename faust_real>
faust_real faust_cu_vec<faust_real>::max()const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val = faust_cu_max(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (faust_real)0.0;
   }
}

template <typename faust_real>
faust_real faust_cu_vec<faust_real>::min()const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val = faust_cu_min(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (faust_real)0.0;
   }
}

template <typename faust_real>
faust_real faust_cu_vec<faust_real>::sum()const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val = faust_cu_sum(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (faust_real)0.0;
   }
}

template <typename faust_real>
faust_real faust_cu_vec<faust_real>::norm() const
{
   if(dim>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val = faust_cu_norm(data, dim);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as vector is empty
      return (faust_real)0.0;
   }
}

template <typename faust_real>
faust_real faust_cu_vec<faust_real>::dot(const faust_cu_vec<faust_real>& cu_v, cublasHandle_t cublasHandle) const
{
//return dot(*this, cu_v, cublasHandle);
   if(size() != cu_v.size())
      handleError("LinAlgebra_cu","dot : the two vectors don't have the same size");
   if(size() > 0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      const faust_cu_vec<faust_real>* cu_v_ptr = &cu_v;
      if(cu_v.getDevice() != device)
         cu_v_ptr = new faust_cu_vec<faust_real>(cu_v, device);

      faust_real result=0.0;

      faust_cu_dot(cublasHandle, size(), getData(), 1, cu_v_ptr->getData(), 1, &result);

      if(cu_v.getDevice() != device)
         delete cu_v_ptr;
      faust_cudaSetDevice(currentGPU);
      return result;
   }
   else
   {
      // display warning as vector is empty
      return (faust_real)0.0;
   }


}

template <typename faust_real>
bool faust_cu_vec<faust_real>::operator==(const faust_cu_vec<faust_real>& cu_v)const
{
   if(data!=NULL && data==cu_v.data && device== cu_v.device)
   {
      if(dim!=cu_v.dim)
         handleError(class_name,"operator== : same data and device but different dimensions");
      return true;     
   }
   else
      return false;
}

template <typename faust_real>
void faust_cu_vec<faust_real>::operator+=(const faust_cu_vec<faust_real>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(class_name,"operator+= : vector dimensions must agree");
   }  
   if(dim == 0)
   {
      // display warning as vector is empty
      return; 
   }   
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   faust_cu_vec vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_add(data, vec_tmp.data, dim);
   }
   else
      kernel_add(data, cu_v.data, dim);

   faust_cudaSetDevice(currentGPU);
}
template <typename faust_real>
void faust_cu_vec<faust_real>::operator-=(const faust_cu_vec<faust_real>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(class_name,"operator-= : vector dimensions must agree");
   }  
   if(dim == 0)
   {
      // display warning as vector is empty
      return; 
   }   
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   faust_cu_vec vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_sub(data, vec_tmp.data, dim);
   }
   else
      kernel_sub(data, cu_v.data, dim);
   faust_cudaSetDevice(currentGPU);
}
template <typename faust_real>
void faust_cu_vec<faust_real>::operator*=(const faust_cu_vec<faust_real>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(class_name,"operator*= : vector dimensions must agree");
   }  
   if(dim == 0)
   {
      // display warning as vector is empty
      return; 
   }   
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);
   faust_cu_vec vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_mult(data, vec_tmp.data, dim);
   }
   else
      kernel_mult(data, cu_v.data, dim);
   faust_cudaSetDevice(currentGPU);
}
template <typename faust_real>
void faust_cu_vec<faust_real>::operator/=(const faust_cu_vec<faust_real>& cu_v)
{
   if(cu_v.dim != dim)
   {
       handleError(class_name,"operator/= : vector dimensions must agree");
   }  
   if(dim == 0)
   {
      // display warning as vector is empty
      return; 
   }   
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   faust_cu_vec vec_tmp;
   if(device != cu_v.device)
   {
      vec_tmp.init(cu_v, device);
      kernel_div(data, vec_tmp.data, dim);
   }
   else
      kernel_div(data, cu_v.data, dim);

   faust_cudaSetDevice(currentGPU);
}


template <typename faust_real>
void faust_cu_vec<faust_real>::operator+=(const faust_real alpha)
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

template <typename faust_real>
void faust_cu_vec<faust_real>::operator-=(const faust_real alpha)
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

template <typename faust_real>
void faust_cu_vec<faust_real>::operator*=(const faust_real alpha)
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

template <typename faust_real>
void faust_cu_vec<faust_real>::operator/=(const faust_real alpha)
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



template <typename faust_real>
void faust_cu_vec<faust_real>::setValues(const faust_real val)
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

template <typename faust_real>
void faust_cu_vec<faust_real>::setValues(const faust_real val, const faust_unsigned_int length)
{
   resize(length);
   setValues(val);
}

template <typename faust_real>
void faust_cu_vec<faust_real>::setZeros()
{setValues(0.0);}

template <typename faust_real>
void faust_cu_vec<faust_real>::setZeros(const faust_unsigned_int length)
{setValues(0.0, length);}

template <typename faust_real>
void faust_cu_vec<faust_real>::setOnes()
{setValues(1.0);}

template <typename faust_real>
void faust_cu_vec<faust_real>::setOnes(const faust_unsigned_int length)
{setValues(1.0, length);}


template <typename faust_real>
faust_real faust_cu_vec<faust_real>::mean_relative_error(const faust_cu_vec<faust_real>& v_ref)const
{
   if(v_ref.size() != size())
     handleError(class_name,"relative_error : sizes are different");

   if(dim == 0)
   {
      // display warning as vector is empty
      return (faust_real)0.0; 
   }   
   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   faust_cu_vec tmp(size());
   kernel_relative_error(tmp.getData(), v_ref.getData(), data, size());

   faust_cudaSetDevice(currentGPU);
   return tmp.mean();
}

template <typename faust_real>
void faust_cu_vec<faust_real>::Display()const
{
   faust_vec<faust_real> v;
   faust_cu2faust(v,*this);
   v.Display();
}

template <typename faust_real>
void faust_cu_vec<faust_real>::print_file(const char* filename)const
{
   faust_vec<faust_real> v;
   faust_cu2faust(v,*this);
   v.print_file(filename);
}

#endif
