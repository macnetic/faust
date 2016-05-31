#ifndef __MatDense_HPP__
#define __MatDense_HPP__

#include "faust_Vect_gpu.h"
#include "faust_MatDense.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifdef __COMPILE_SPMAT__
	#include "faust_MatSparse.h"
	#include "faust_MatSparse_gpu.h"
#endif
#include "faust_reduce_gpu.h"
#include "kernels.h"
#include "faust_cuda.h"
#include "faust_gpu2cpu.h"

#include "linear_algebra_gpu.h"


using namespace std;

template <typename FPP>
const char * Faust::MatDense<FPP,Gpu>::class_name = "Faust::MatDense<FPP,Gpu>::";

template <typename FPP>
Faust::MatDense<FPP,Gpu>::MatDense(const FPP *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.start();
   else
      t_constructor_from_host.start();
#endif

   if(nbRow*nbCol == 0)
   {
      resize(nbRow, nbCol, dstDevice);
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.stop();
   else
      t_constructor_from_host.stop();
#endif
      return;
    }
    if(data_ == NULL)
            handleError(class_name, "Faust::MatDense(const faust_real*,const faust_unsigned_int, const faust_unsigned_int, bool, int=FAUST_DEFAULT_CUDA_DEVICE, int=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t=0) : data pointer argument is NULL");

    if(dataFromGPU)
        copyFromDevice(data_, nbRow, nbCol, dstDevice, srcDevice, stream);
    else
        copyFromHost(data_, nbRow, nbCol, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.stop();
   else
      t_constructor_from_host.stop();
#endif
}


//WARNING : if cu_M is Identity, the programm will failed (so an error is launched)
template <typename FPP>
Faust::MatDense<FPP,Gpu>::MatDense(const Faust::MatDense<FPP,Gpu>& cu_M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.start();
#endif
    //update the flag before resizing the matrix
    if (cu_M.isIdentity || cu_M.isZeros)
	handleError(class_name, "Faust::MatDense(const Faust::MatDense<..,Gpu>& cu_M, ...) :  cu_M is Identity, Zeros are not handle by the constructor");
    if(cu_M.isIdentity || cu_M.isZeros || cu_M.dim1*cu_M.dim2==0)
    {
	resize(cu_M.dim1, cu_M.dim2, dstDevice);
	if(cu_M.isIdentity)
            isIdentity = true;
        if(cu_M.isZeros)
            isZeros = true;

    }
    else
        copyFromDevice(cu_M.data, cu_M.dim1, cu_M.dim2, dstDevice, cu_M.device, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.stop();
#endif
}

template <typename FPP>
Faust::MatDense<FPP,Gpu>::MatDense(const Faust::MatDense<FPP,Cpu>& M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL),device(FAUST_DEFAULT_CUDA_DEVICE)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.start();
#endif
    if(M.estIdentite() || M.estNulle() || M.getNbRow()*M.getNbCol()==0)
    {
        resize(M.getNbRow(), M.getNbCol(), dstDevice);
   dim1 = M.getNbRow();
        dim2 = M.getNbCol();
        if(M.estIdentite())
            isIdentity = true;
        if(M.estNulle())
        isZeros = true;
    }
    else
        copyFromHost(M.getData(), M.getNbRow(), M.getNbCol(), dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.stop();
#endif
}

#ifdef __COMPILE_SPMAT__
template <typename FPP>
Faust::MatDense<FPP,Gpu>::MatDense(const Faust::MatSparse<FPP,Gpu>& cu_S,Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.start();
#endif
    if(cu_S.getNonZeros() == 0)
    {
        if(cu_S.getRowPtr()!=NULL || cu_S.getColInd()!=NULL || cu_S.getValues()!=NULL)
            handleError(class_name, "Faust::MatDense(const faust_cu_spmat&, ...) : rowPtr,  colInd or values pointer is not NULL");

        setZeros(cu_S.getNbRow(), cu_S.getNbCol(), dstDevice);
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.stop();
#endif
        return;
    }

    resize(cu_S.getNbRow(), cu_S.getNbCol(), dstDevice);

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    const Faust::MatSparse<FPP,Gpu>* cu_S_dst = &cu_S;
    if(dstDevice != cu_S.getDevice())
        cu_S_dst = new Faust::MatSparse<FPP,Gpu>(cu_S, dstDevice);

    faust_cu_csr2dense(spblasHandle.GetHandle(),
            dim1,dim2,
            cu_S_dst->getDescr(),
            cu_S_dst->getValues(), cu_S_dst->getRowPtr(), cu_S_dst->getColInd(),
            data, dim1);


    if(dstDevice != cu_S.getDevice())
        delete    cu_S_dst;
    cu_S_dst = NULL;

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.stop();
#endif
}
#endif

template <typename FPP>
Faust::MatDense<FPP,Gpu>::MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_size.start();
#endif
    resize(nbRow, nbCol, dstDevice);
#ifdef __COMPILE_TIMERS__
t_constructor_from_size.stop();
#endif
}



/// GETTEUR SETTEUR ///

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::_create(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_)
{
#ifdef __COMPILE_TIMERS__
t_create.start();
#endif
   if(nbRow<0 || nbCol<0)
       handleError(class_name, "_create : incorrect dimensions");


   if(nbRow*nbCol==0)
      data = NULL;

   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device_);

      if (data == NULL)
         faust_cudaMalloc((void**)&data, (nbRow*nbCol)*sizeof(FPP));
      else
         handleError(class_name, "_create : data has already been allocated on GPU");


      faust_cudaSetDevice(currentGPU);
   }
   dim1 = nbRow;
   dim2 = nbCol;
   device = device_;
   isZeros = false;
   isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_create.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::_clear()
{
#ifdef __COMPILE_TIMERS__
t_clear.start();
#endif
   if(dim1*dim2>0 && !isZeros && !isIdentity)
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
      handleError(class_name, "_clear : data of empty matrix, identity matrix or zeros matrix should be NULL");

   dim1 = 0;
   dim2 = 0;
   isIdentity = false;
   isZeros = false;
   data = NULL;

#ifdef __COMPILE_TIMERS__
t_clear.stop();
#endif
}


template <typename FPP>
void Faust::MatDense<FPP,Gpu>::resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_)
{
   if(nbRow*nbCol!=dim1*dim2 || (device_!=device) || isIdentity || isZeros)
   {
      _clear();
      _create(nbRow, nbCol, device_);
   }
   else
   {
      dim1 = nbRow;
      dim2 = nbCol;
      isIdentity = false;
      isZeros = false;
      device = device_;
   }
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::copyFromHost(const FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_copy_from_host.start();
#endif
    resize(nbRow, nbCol, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(data, data_, nbRow*nbCol*sizeof(FPP), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_copy_from_host.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::copyFromDevice(const FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_copy_from_device.start();
#endif
    resize(nbRow, nbCol, dstDevice);
    faust_cudaMemcpyPeerAsync(data, dstDevice, data_, srcDevice, nbRow*nbCol*sizeof(FPP), stream);
#ifdef __COMPILE_TIMERS__
t_copy_from_device.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::copyToHost(FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream/*=0*/)const
{
#ifdef __COMPILE_TIMERS__
t_copy_to_host.start();
#endif
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(data_, data, nbRow*nbCol*sizeof(FPP), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_copy_to_host.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::copyToDevice(FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
#ifdef __COMPILE_TIMERS__
t_copy_to_device.start();
#endif
       faust_cudaMemcpyPeerAsync(data_, dstDevice, data, device, nbRow*nbCol*sizeof(FPP), stream);
#ifdef __COMPILE_TIMERS__
t_copy_to_device.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::init(const Faust::MatDense<FPP,Gpu>& cu_M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_init.start();
#endif
    if(cu_M.isIdentity || cu_M.isZeros || data==NULL )
    {
        resize(0,0);
        dim1 = cu_M.dim1;
        dim2 = cu_M.dim2;
        isIdentity = cu_M.isIdentity;
        isZeros = cu_M.isZeros;
        data = NULL ;
        device = dstDevice;
    }
    else
        copyFromDevice(cu_M.data, cu_M.dim1, cu_M.dim2, dstDevice, cu_M.device, stream);
#ifdef __COMPILE_TIMERS__
t_init.stop();
#endif

}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_move_to_device.start();
#endif
   if(device == dstDevice)
   {

#ifdef __COMPILE_TIMERS__
t_move_to_device.stop();
#endif
      return;
   }
   if (data == NULL)
   {
      device = dstDevice;
#ifdef __COMPILE_TIMERS__
t_move_to_device.stop();
#endif
      return;
   }

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(dstDevice);

   FPP* data_tmp;
   faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(FPP));
   faust_cudaMemcpyPeerAsync(data_tmp, dstDevice, data, device, dim1*dim2*sizeof(FPP), stream);
   faust_cudaFree(data);
   data = data_tmp;
   device = dstDevice;

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_move_to_device.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_)
{
#ifdef __COMPILE_TIMERS__
t_set_zeros.start();
#endif
    resize(0,0,device_);
    dim1 = nbRow;
    dim2 = nbCol;
    isZeros = true;
    isIdentity = false;
#ifdef __COMPILE_TIMERS__
t_set_zeros.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::setEyes(const faust_unsigned_int nbRow, const int device_)
{
#ifdef __COMPILE_TIMERS__
t_set_eyes.start();
#endif
    if(nbRow==device_)
        cerr << __FILE__ << ":" << __LINE__ << " : Warning - prototype is setEyes(const faust_unsigned_int nbRow, const int device_) and not setEyes(const int nbRow, const int nbCol). Please check whether syntax is correct" << endl;
    setZeros(nbRow, nbRow, device_);
    isIdentity = true;
    isZeros = false;
#ifdef __COMPILE_TIMERS__
t_set_eyes.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::setEyes()
{
   if(dim1!=dim2)
      handleError(class_name,"setEyes() : GPU matrix is not square");
   setEyes(dim1, device);
}
 /// OPERATION BASIQUE ///

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::transpose(Faust::BlasHandle<Gpu> blasHandle)
{

#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif

   if(isZeros)
   {
      int dim1_old = dim1;
      int dim2_old = dim2;
      resize(0,0);
      dim1 = dim2_old;
      dim2 = dim1_old;
      isZeros = true;
      #ifdef __COMPILE_TIMERS__
         t_transpose.stop();
      #endif
      return;
   }
   if(isIdentity)
   {
      #ifdef __COMPILE_TIMERS__
         t_transpose.stop();
      #endif
      return;
   }
   if (dim1*dim2 > 0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      FPP* data_tmp;
      faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(FPP));
      faust_cudaMemcpy(data_tmp, data, dim1*dim2*sizeof(FPP), cudaMemcpyDeviceToDevice);

      FPP alpha=1.0f;
      FPP beta=0.0f;
      // According to the cublas doc, in-place mode if C = B, ldc = ldb and transb = CUBLAS_OP_N
      faust_cu_geam(blasHandle.GetHandle(),
         CUBLAS_OP_T, CUBLAS_OP_N,
         dim2, dim1,
         &alpha,
         data, dim1,
         &beta,
         data_tmp, dim2,
         data_tmp, dim2);


      faust_cudaFree(data);
      data = data_tmp;

      faust_cudaSetDevice(currentGPU);
   }
   int dim_tmp = dim1;
   dim1 = dim2;
   dim2 = dim_tmp;

#ifdef __COMPILE_TIMERS__
t_transpose.stop();
#endif

}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::init_from_transpose(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::BlasHandle<Gpu> blasHandle)
{
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.start();
#endif

if(cu_A == (*this))
            handleError(class_name, "init_from_transpose(const Faust::MatDense<FPP,Gpu>&) : input GPU data is the same that the current GPU data. Try using Faust::MatDense<FPP,Gpu>::transpose instead, in this particular case.");

   if(cu_A.isZeros)
   {
      setZeros(cu_A.dim2, cu_A.dim1);
      #ifdef __COMPILE_TIMERS__
         t_init_from_transpose.stop();
      #endif
      return;
   }
   if(cu_A.isIdentity)
   {
      setEyes(cu_A.dim1);
      #ifdef __COMPILE_TIMERS__
         t_init_from_transpose.stop();
      #endif
      return;
   }
   resize(cu_A.dim2, cu_A.dim1);
   if (cu_A.dim1*cu_A.dim2 > 0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      FPP* data_tmp = NULL;
      if(device == cu_A.device)
          data_tmp = cu_A.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, cu_A.dim1*cu_A.dim2*sizeof(FPP));
         faust_cudaMemcpyPeer(data_tmp, device, cu_A.data, cu_A.device, cu_A.dim1*cu_A.dim2*sizeof(FPP));
      }

      FPP alpha=1.0f;
      FPP beta=0.0f;
      faust_cu_geam(blasHandle.GetHandle(),
         CUBLAS_OP_T, CUBLAS_OP_N,
         cu_A.dim2, cu_A.dim1,
         &alpha,
         data_tmp, cu_A.dim1,
         &beta,
         data, cu_A.dim2,
         data, cu_A.dim2);


      if(device != cu_A.device)
      {
         faust_cudaFree(data_tmp);
         data_tmp = NULL;
      }


      faust_cudaSetDevice(currentGPU);
   }

#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif

}

#ifdef __COMPILE_SPMAT__
template <typename FPP>
void Faust::MatDense<FPP,Gpu>::init_from_transpose(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.start();
#endif

    if(cu_S.getNbRow()*cu_S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.getNbRow(), cu_S.getNbCol());
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif
       return;
    }

    if(cu_S.getNonZeros() == 0)
    {
       // display warning as matrix is zeros
       setZeros(cu_S.getNbRow(), cu_S.getNbCol());
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif
       return;
    }


   resize(cu_S.getNbCol(), cu_S.getNbRow());
   if (cu_S.getNbRow()*cu_S.getNbCol() > 0)
   {
      if (cu_S.getNonZeros()==0)
         handleError(class_name, "init_from_transpose(faust_cu_spmat, ...) : impossible to create full matrix from uninitialized sparse matrix");
      else
      {
         int currentGPU;
         faust_cudaGetDevice(&currentGPU);
         faust_cudaSetDevice(device);

         if(device == cu_S.getDevice())
            faust_cu_csc2dense(spblasHandle.GetHandle(),
               cu_S.getNbCol(), cu_S.getNbRow(), cu_S.getDescr(),
               cu_S.getValues(), cu_S.getColInd(), cu_S.getRowPtr(),
               data, cu_S.getNbCol());
         else
         {
            int *rowptr_tmp, *colind_tmp;
            FPP *values_tmp;
            faust_cudaMalloc((void**)&rowptr_tmp, (cu_S.getNbRow()+1)*sizeof(int));
            faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
            faust_cudaMalloc((void**)&values_tmp, cu_S.getNonZeros()*sizeof(FPP));
            faust_cudaMemcpyPeer(rowptr_tmp, device, cu_S.getRowPtr(), cu_S.getDevice(), (cu_S.getNbRow()+1)*sizeof(int));
            faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
            faust_cudaMemcpyPeer(values_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(FPP));
            faust_cu_csc2dense(spblasHandle.GetHandle(),
               cu_S.getNbCol(), cu_S.getNbRow(), cu_S.getDescr(),
               values_tmp, colind_tmp, rowptr_tmp,
               data, cu_S.getNbCol());
            faust_cudaFree(rowptr_tmp);
            faust_cudaFree(colind_tmp);
            faust_cudaFree(values_tmp);
         }

          faust_cudaSetDevice(currentGPU);
      }
   }
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif

}
#endif



template <typename FPP>
FPP Faust::MatDense<FPP,Gpu>::max() const
{
#ifdef __COMPILE_TIMERS__
t_max.start();
#endif
   if(isZeros)
   {
#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif
      return (FPP)0.0;
   }
   else if (isIdentity)
   {
#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif
      return (FPP)1.0;
   }
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val =  faust_cu_max(data, dim1*dim2);
      faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif
      return val;
   }
   else
   {
      // display warning as matrix is empty
#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif
      return (FPP)0.0;
   }
#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif
}

template <typename FPP>
FPP Faust::MatDense<FPP,Gpu>::min() const
{
#ifdef __COMPILE_TIMERS__
t_min.start();
#endif
   if(isZeros || isIdentity)
   {
#ifdef __COMPILE_TIMERS__
t_min.stop();
#endif
      return (FPP)0.0;
   }
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val =  faust_cu_min(data, dim1*dim2);
      faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_min.stop();
#endif
      return val;
   }
   else
   {
      // display warning as matrix is empty
#ifdef __COMPILE_TIMERS__
t_min.stop();
#endif
      return (FPP)0.0;
   }
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::abs()
{
#ifdef __COMPILE_TIMERS__
t_abs.start();
#endif
   if(isZeros || isIdentity)
   {
#ifdef __COMPILE_TIMERS__
t_abs.stop();
#endif
      return ;
   }
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      kernel_abs(data, dim1*dim2);
      faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_abs.stop();
#endif
      return;
   }
   else
   {
#ifdef __COMPILE_TIMERS__
t_abs.stop();
#endif
      // display warning as matrix is empty
      return;
   }
#ifdef __COMPILE_TIMERS__
t_abs.stop();
#endif
}

template <typename FPP>
FPP Faust::MatDense<FPP,Gpu>::norm() const
{
#ifdef __COMPILE_TIMERS__
t_norm.start();
#endif
   if(isZeros)
   {
#ifdef __COMPILE_TIMERS__
t_norm.stop();
#endif
      return (FPP)0.0;
   }
   else if(isIdentity)
   {
#ifdef __COMPILE_TIMERS__
t_norm.stop();
#endif
      return (FPP)sqrt(dim1);
   }
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val =  faust_cu_norm(data, dim1*dim2);
      faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_norm.stop();
#endif
      return val;
   }
   else
   {
      // display warning as matrix is empty
#ifdef __COMPILE_TIMERS__
t_norm.stop();
#endif
      return (FPP)0.0;
   }
#ifdef __COMPILE_TIMERS__
t_norm.stop();
#endif
}

template <typename FPP>
FPP Faust::MatDense<FPP,Gpu>::trace() const
{
#ifdef __COMPILE_TIMERS__
t_trace.start();
#endif
   if(dim1 != dim2)
      handleError(class_name, "norm : matrix must be square");

   if (isZeros)
   {
#ifdef __COMPILE_TIMERS__
t_trace.stop();
#endif
      return (FPP)0.0;
   }
   else if (isIdentity)
   {
#ifdef __COMPILE_TIMERS__
t_trace.stop();
#endif
      return (FPP)dim1;
   }
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      Faust::Vect<FPP,Gpu> cu_diag(dim1);
      kernel_get_diag(cu_diag.getData(), data, dim1);
      faust_cudaSetDevice(currentGPU);

#ifdef __COMPILE_TIMERS__
t_trace.stop();
#endif
      return faust_cu_sum(cu_diag.getData(), dim1);
   }
   else
   {
      // display warning as matrix is empty
#ifdef __COMPILE_TIMERS__
t_trace.stop();
#endif
      return (FPP)0.0;
   }
#ifdef __COMPILE_TIMERS__
t_trace.stop();
#endif
}





template <typename FPP>
 FPP Faust::MatDense<FPP,Gpu>::spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag, Faust::BlasHandle<Gpu> blasHandle) const
{
   #ifdef __COMPILE_TIMERS__
      t_spectral_norm2.start();
   #endif
   if(isZeros)
   {
      flag = -2;
      #ifdef __COMPILE_TIMERS__
      t_spectral_norm2.stop();
      #endif
      return 0;
   }

   if(isIdentity)
   {
      flag = -3;
      #ifdef __COMPILE_TIMERS__
      t_spectral_norm2.stop();
      #endif
      return 1;
   }

   faust_unsigned_int nb_row = getNbRow();
   faust_unsigned_int nb_col = getNbCol();


   Faust::MatDense<FPP,Gpu> AtA;
   if (nb_row <= nb_col)
      gemm<FPP>((*this),(*this),AtA,1.,0,'N','T',blasHandle);
   else
      gemm<FPP>((*this),(*this),AtA,1.,0,'T','N', blasHandle);



   FPP  res=std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag,blasHandle));



   #ifdef __COMPILE_TIMERS__
      t_spectral_norm2.stop();
   #endif
   return res;

}


template <typename FPP>
void Faust::MatDense<FPP,Gpu>::scalarMultiply(const FPP lambda)
{
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.start();
#endif
   if(dim1*dim2 == 0)
   {
      // display warning as matrix is empty
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
      return;
   }

   if(isZeros)
   {
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      if(isIdentity)
      {
         // display warning as it could be better to create an instance of faust_cu_spmat class instead

         resize(dim1,dim2);
         kernel_memset(data, (FPP)0.0, dim1*dim2);
         kernel_add_diag_const(data, lambda, dim1);
      }
      else
      {
         kernel_mult_const(data, lambda, dim1*dim2);
      }

      faust_cudaSetDevice(currentGPU);
   }

#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
}

template <typename FPP>
bool Faust::MatDense<FPP,Gpu>::operator==(const Faust::MatDense<FPP,Gpu>& cu_M)const
{
   if(data!=NULL && data==cu_M.data && device==cu_M.device)
   {
      if(dim1!=cu_M.dim1 || dim2!=cu_M.dim2 && !isIdentity && !isZeros)
         handleError(class_name,"operator== : same data and device but different dimensions or identity or zeros");
      return true;
   }
   else
      return false;
}

    /// SURCHARGE OPERATEUR ///
  // affectation
template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator=(const Faust::MatDense<FPP,Gpu>& cu_M)
{
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.start();
#endif
   if(cu_M.dim1*cu_M.dim2 == 0)
   {
      // display warning as matrix is empty
      resize(cu_M.dim1, cu_M.dim2);
      data = NULL ;
      device = cu_M.device;
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.stop();
#endif
      return;
   }
   if(cu_M.isIdentity || cu_M.isZeros)
   {
      resize(0,0);
      dim1 = cu_M.dim1;
      dim2 = cu_M.dim2;
      isIdentity = cu_M.isIdentity;
      isZeros = cu_M.isZeros;
      data = NULL ;
      device = cu_M.device;
   }
   else
   {
       copyFromDevice(cu_M.data, cu_M.dim1, cu_M.dim2, cu_M.device, cu_M.device);
   }
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator=(const Faust::MatDense<FPP,Cpu>& M)
{
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.start();
#endif

    if(M.getNbRow()*M.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(M.getNbRow(), M.getNbCol());
       data = NULL ;
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.stop();
#endif
       return;
    }
    else if(M.estIdentite() || M.estNulle())
    {
        resize(0,0);
        dim1 = M.getNbRow();
        dim2 = M.getNbCol();
        isIdentity = M.estIdentite();
        isZeros = M.estNulle();
        data = NULL ;
    }
    else
        copyFromHost(M.getData(), M.getNbRow(), M.getNbCol(), device);
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.stop();
#endif
}

#ifdef __COMPILE_SPMAT__
template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator=(const Faust::MatSparse<FPP,Cpu>& S)
{
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.start();
#endif
    if(S.getNbRow()*S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(S.getNbRow(), S.getNbCol());
       data = NULL ;
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.stop();
#endif
       return;
    }
    if(S.getNonZeros()==0)
    {
        setZeros(S.getNbRow(), S.getNbCol());
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.stop();
#endif
        return;
    }

    Faust::MatDense<FPP,Cpu> M;
    M=S;

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpy(data, M.getData(), dim1*dim2*sizeof(FPP), cudaMemcpyHostToDevice);

    faust_cudaSetDevice(currentGPU);

    isZeros = false;
    isIdentity = false;
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_host.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator=(const Faust::MatSparse<FPP,Gpu>& S)
{handleError(class_name, "Faust::MatDense<FPP,Gpu>::operator=(const faust_cu_spmat&) is not defined. Use Faust::MatDense<FPP,Gpu>::init_from_cu_spmat(const faust_cu_spmat&, cusparseHandle_t) instead");}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::init_from_cu_spmat(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle, const FPP coeff /*=1.0*/)
{
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.start();
#endif
    if(cu_S.getNbRow()*cu_S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.getNbRow(), cu_S.getNbCol());
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.stop();
#endif
       return;
    }

    if(cu_S.getNonZeros() == 0)
    {
       // display warning as matrix is zeros
       setZeros(cu_S.getNbRow(), cu_S.getNbCol());
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.stop();
#endif
       return;
    }


    resize(cu_S.getNbRow(), cu_S.getNbCol());

   if(cu_S.getNbRow()*cu_S.getNbCol()>0)
   {
         int currentGPU;
         faust_cudaGetDevice(&currentGPU);
         faust_cudaSetDevice(device);

         if(device==cu_S.getDevice() && coeff==1.0 )
         {
               faust_cu_csr2dense(spblasHandle.GetHandle(),
                  cu_S.getNbRow(), cu_S.getNbCol(), cu_S.getDescr(),
                  cu_S.getValues(), cu_S.getRowPtr(), cu_S.getColInd(),
                  data, cu_S.getNbRow());
         }
         else
         {
            FPP *values_tmp;
            faust_cudaMalloc((void**)&values_tmp, cu_S.getNonZeros()*sizeof(FPP));
            faust_cudaMemcpyPeer(values_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(FPP));
            if(coeff != 1.0)
               kernel_mult_const(values_tmp, coeff, cu_S.getNonZeros());

            if(device!=cu_S.getDevice())
            {
               int *rowptr_tmp, *colind_tmp;

               faust_cudaMalloc((void**)&rowptr_tmp, cu_S.getNonZeros()*sizeof(int));
               faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));

               faust_cudaMemcpyPeer(rowptr_tmp, device, cu_S.getRowPtr(), cu_S.getDevice(), (cu_S.getNbRow()+1)*sizeof(int));
               faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));


               faust_cu_csr2dense(spblasHandle.GetHandle(),
                     cu_S.getNbRow(), cu_S.getNbCol(), cu_S.getDescr(),
                     values_tmp, rowptr_tmp, colind_tmp,
                     data, cu_S.getNbRow());

               faust_cudaFree(rowptr_tmp);
               faust_cudaFree(colind_tmp);
               rowptr_tmp = NULL;
               colind_tmp = NULL;
            }
            else
            {
               faust_cu_csr2dense(spblasHandle.GetHandle(),
                     cu_S.getNbRow(), cu_S.getNbCol(), cu_S.getDescr(),
                     values_tmp, cu_S.getRowPtr(), cu_S.getColInd(),
                     data, cu_S.getNbRow());

            }
            faust_cudaFree(values_tmp);
            values_tmp = NULL;
         }

         faust_cudaSetDevice(currentGPU);


   }

/*----------- ancienne version (sans appel a kernel_sparse2full) -------------------------------
    faust_unsigned_int* rowind = new faust_unsigned_int[cu_S.getNonZeros()];
    faust_unsigned_int* colind = new faust_unsigned_int[cu_S.getNonZeros()];
    faust_real* values = new faust_real[cu_S.getNonZeros()];
    faust_real* ptr_data = new faust_real[dim1*dim2];
    for (int i=0 ; i<cu_S.getNonZeros() ; i++)
        rowind[i]=0;
    for (int i=0 ; i<cu_S.getNonZeros() ; i++)
        colind[i]=0;
    for (int i=0 ; i<cu_S.getNonZeros() ; i++)
        values[i]=0.0;
    for (int i=0 ; i<dim1*dim2 ; i++)
        ptr_data[i]=0.0;

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(cu_S.getDevice());
    faust_cudaMemcpy(rowind, cu_S.getRowind(), cu_S.getNonZeros*sizeof(faust_unsigned_int), cudaMemcpyDeviceToHost);
    faust_cudaMemcpy(colind, cu_S.getColInd(), cu_S.getNonZeros*sizeof(faust_unsigned_int), cudaMemcpyDeviceToHost);
    faust_cudaMemcpy(values, cu_S.getValues(), cu_S.getNonZeros*sizeof(faust_real), cudaMemcpyDeviceToHost);

    if(cu_S.getIdxBase() == CUSPARSE_INDEX_BASE_ZERO)
        for(int i=0 ; i< cu_S.getNonZeros() ; i++)
            ptr_data[colind[i] * dim1 + rowind[i]] = values[i];
    else (cu_S.getIdxBase() == CUSPARSE_INDEX_BASE_ONE)
        for(int i=0 ; i< cu_S.getNonZeros() ; i++)
            ptr_data[(colind[i]-1) * dim1 + (rowind[i]-1)] = values[i];

    faust_cudaSetDevice(device);
    faust_cudaMemcpy(data, ptr_data, dim1*dim2*sizeof(faust_real), cudaMemcpyHostToDevice);
    faust_cudaSetDevice(currentGPU);

    delete[] rowind ; rowind=NULL;
    delete[] colind ; colind=NULL;
    delete[] values ; values=NULL;
    delete[] ptr_data ; ptr_data=NULL;
------------------------------------------------------------------------------------------------*/
#ifdef __COMPILE_TIMERS__
t_operator_equal_from_device.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator+=(const Faust::MatSparse<FPP,Gpu>& cu_S)
{handleError(class_name, "Faust::MatDense<FPP,Gpu>::operator+=(const faust_cu_spmat&) is not defined. Use Faust::MatDense<FPP,Gpu>::add(const faust_cu_spmat&, Faust::SpBlasHandle<Gpu>) instead");}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::add(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
t_add_cuspmat.start();
#endif
   if(!(isZeros || isIdentity || data))
      handleError(class_name,"add : uninitialized matrix");

   if(dim1!=cu_S.getNbRow() || dim2!=cu_S.getNbCol())
      handleError(class_name,"add : incorrect matrix dimensions");

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   if(isZeros && cu_S.getNonZeros()==0)
      setZeros(dim1, dim2);
   else if (isZeros)
      init_from_cu_spmat(cu_S, spblasHandle);
   else if (cu_S.getNonZeros()==0){} // We don't do anything but it is written not to go to else
   else if(isIdentity)
   {
      init_from_cu_spmat(cu_S, spblasHandle);
      kernel_add_diag_const(data, (FPP)1.0, dim1);
   }
   else
   {

      faust_cudaSetDevice(cu_S.getDevice());

      int* rowind_tmp_src;
      faust_cudaMalloc((void**)&rowind_tmp_src, cu_S.getNonZeros()*sizeof(int));
      faust_cu_csr2coo(spblasHandle.GetHandle(), cu_S.getRowPtr(),
         cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaSetDevice(device);


      if(device == cu_S.getDevice())
      {
         faust_cu_csr2coo(spblasHandle.GetHandle(), cu_S.getRowPtr(),
            cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

         kernel_add_sparse2full(data, rowind_tmp_src, cu_S.getColInd(), cu_S.getValues(), cu_S.getNonZeros(), dim1);

      }
      else
      {
         int *rowind_tmp_dst, *colind_tmp;
         FPP* data_tmp;

         faust_cudaMalloc((void**)&rowind_tmp_dst, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&data_tmp, cu_S.getNonZeros()*sizeof(FPP));
         faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(rowind_tmp_dst, device, rowind_tmp_src, cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(data_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(FPP));
         kernel_add_sparse2full(data, rowind_tmp_dst, colind_tmp, data_tmp, cu_S.getNonZeros(), dim1);
         faust_cudaFree(rowind_tmp_dst);
         faust_cudaFree(colind_tmp);
         faust_cudaFree(data_tmp);
      }
      faust_cudaFree(rowind_tmp_src);
   }

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_add_cuspmat.stop();
#endif
}


template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator-=(const Faust::MatSparse<FPP,Gpu>& cu_S)
{
   handleError(class_name, "Faust::MatDense<FPP,Gpu>::operator-=(const faust_cu_spmat&) is not defined. Use Faust::MatDense<FPP,Gpu>::sub(const faust_cu_spmat&, cusparseHandle_t) instead");
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::sub(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
t_sub.start();
#endif
   if(!(isZeros || isIdentity || data))
      handleError(class_name,"sub : uninitialized matrix");

   if(dim1!=cu_S.getNbRow() || dim2!=cu_S.getNbCol())
      handleError(class_name,"sub : incorrect matrix dimensions");

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   if(isZeros && cu_S.getNonZeros()==0)
      setZeros(dim1, dim2);
   else if (isZeros)
      init_from_cu_spmat(cu_S, spblasHandle, -1.0);
   else if (cu_S.getNonZeros()==0){} // We don't do anything but it is written not to go to else
   else if(isIdentity)
   {
      init_from_cu_spmat(cu_S,spblasHandle, -1.0);
      kernel_add_diag_const(data, (FPP)1.0, dim1);
   }
   else
   {
       // if *this and cu_S are the same object, we can use kernel_add, so there is no need to know if these matrix are the same objet or not

      faust_cudaSetDevice(cu_S.getDevice());

      int* rowind_tmp_src;
      cudaMalloc((void**)&rowind_tmp_src, cu_S.getNonZeros()*sizeof(int));
      faust_cu_csr2coo(spblasHandle.GetHandle(), cu_S.getRowPtr(),
         cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaSetDevice(device);


      if(device == cu_S.getDevice())
      {
         faust_cu_csr2coo(spblasHandle.GetHandle(), cu_S.getRowPtr(),
            cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

         kernel_sub_sparse2full(data, rowind_tmp_src, cu_S.getColInd(), cu_S.getValues(), cu_S.getNonZeros(), dim1);

      }
      else
      {
         int *rowind_tmp_dst, *colind_tmp;
         FPP* data_tmp;

         faust_cudaMalloc((void**)&rowind_tmp_dst, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&data_tmp, cu_S.getNonZeros()*sizeof(FPP));
         faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(rowind_tmp_dst, device, rowind_tmp_src, cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(data_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(FPP));
         kernel_sub_sparse2full(data, rowind_tmp_dst, colind_tmp, data_tmp, cu_S.getNonZeros(), dim1);
         faust_cudaFree(rowind_tmp_dst);
         faust_cudaFree(colind_tmp);
         faust_cudaFree(data_tmp);
      }
      faust_cudaFree(rowind_tmp_src);
   }

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_sub.stop();
#endif
}
#endif //#ifdef __COMPILE_SPMAT__


template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator+=(const Faust::MatDense<FPP,Gpu>& cu_M)
{
#ifdef __COMPILE_TIMERS__
t_operator_plus_equal.start();
#endif
   if(!((isZeros || isIdentity || data) && (cu_M.isZeros || cu_M.isIdentity || cu_M.data)))
      handleError(class_name,"operator+= : uninitialized matrix");

   if(dim1!=cu_M.dim1 || dim2!=cu_M.dim2)
      handleError(class_name,"operator+= : incorrect matrix dimensions");

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   if(isZeros && cu_M.isZeros)
      setZeros(dim1, dim2);
   else if (isZeros)
      this->operator=(cu_M);
   else if (cu_M.isZeros){} // We don't do anything but it is written not to go to else

   else if(isIdentity && cu_M.isIdentity)
   {
      resize(dim1, dim2);
      kernel_memset(data, (FPP)0.0, dim1*dim2);
      kernel_add_diag_const(data, (FPP)2.0, dim1);
   }
   else if(isIdentity)
   {
      this->operator=(cu_M);
      kernel_add_diag_const(data, (FPP)1.0, dim1);
   }
   else if(cu_M.isIdentity)
      kernel_add_diag_const(data, (FPP)1.0, dim1);

   else
   {
       // if *this and cu_M are the same object, we can use kernel_add, so there is no need to know if these matrix are the same objet or not
      if(device == cu_M.device)
         kernel_add(data, cu_M.data, dim1*dim2);
      else
      {
         FPP* data_tmp;
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(FPP));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(FPP));
         kernel_add(data, data_tmp, dim1*dim2);
         faust_cudaFree(data_tmp);
      }
   }

   faust_cudaSetDevice(currentGPU);

#ifdef __COMPILE_TIMERS__
t_operator_plus_equal.stop();
#endif
}



template <typename FPP>
void Faust::MatDense<FPP,Gpu>::operator-=(const Faust::MatDense<FPP,Gpu>& cu_M)
{
#ifdef __COMPILE_TIMERS__
t_operator_less_equal.start();
#endif
   if(!isZeros && !isIdentity && data==NULL)
      handleError(class_name,"operator-= : uninitialized matrix");

   if(dim1!=cu_M.dim1 || dim2!=cu_M.dim2)
      handleError(class_name,"operator-= : incorrect matrix dimensions");

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);



   if(isZeros && cu_M.isZeros)
      setZeros(dim1, dim2);
   else if (cu_M.isZeros){} // We don't do anything but it is written not to go to else
   else if(isIdentity && cu_M.isIdentity)
      setZeros(dim1,dim2);
   else if(cu_M.isIdentity)
      kernel_add_diag_const(data, (FPP)-1.0, dim1);
   else
   {
      FPP* data_tmp = NULL;
      if(device == cu_M.device)
         data_tmp = cu_M.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(FPP));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(FPP));
      }

      if (isZeros)
      {
         resize(dim1, dim2);
         kernel_memset(data, (FPP)0.0, dim1*dim2);
         kernel_sub(data, data_tmp, dim1*dim2);
      }
      else if(isIdentity)
      {
         resize(dim1, dim2);
         kernel_memset(data, (FPP)0.0, dim1*dim2);
         kernel_sub(data, data_tmp, dim1*dim2);
         kernel_add_diag_const(data, (FPP)1.0, dim1);
      }
      else if(device==cu_M.device && data==cu_M.data)
      {
         // if *this and cu_M are the same object, we could use kernel_sub, but in order to use less memory, we set the matrix to zero.
         setZeros(dim1,dim2);
      }
      else
      {
         kernel_sub(data, data_tmp, dim1*dim2);
      }

      if(device != cu_M.device)
      {
         faust_cudaFree(data_tmp);
         data_tmp = NULL;
      }
   }

   faust_cudaSetDevice(currentGPU);

#ifdef __COMPILE_TIMERS__
t_operator_less_equal.stop();
#endif
}



template <typename FPP>
void Faust::MatDense<FPP,Gpu>::scalarMultiply(const Faust::MatDense<FPP,Gpu>& cu_M)
{
#ifdef __COMPILE_TIMERS__
t_hadamard_product.start();
#endif
   if(!isZeros && !isIdentity && data==NULL)
      handleError(class_name,"scalarMultiply : uninitialized matrix");

   if(dim1!=cu_M.dim1 || dim2!=cu_M.dim2)
      handleError(class_name,"scalarMultiply : incorrect matrix dimensions");

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   if(isZeros || cu_M.isZeros)
      setZeros(dim1, dim2);

   else if(isIdentity && cu_M.isIdentity){} // We don't do anything but it is written not to go to else

   else if(cu_M.isIdentity)
   {
      Faust::MatDense<FPP,Gpu> cu_A_tmp(*this);
      kernel_memset(data, (FPP)0.0, dim1*dim2);
      kernel_copy_diag(data, cu_A_tmp.data, dim1);
   }
   else
   {
      FPP* data_tmp = NULL;
      if(device == cu_M.device)
         data_tmp = cu_M.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(FPP));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(FPP));
      }

      if(isIdentity)
      {
         kernel_memset(data, (FPP)0.0, dim1*dim2);
         kernel_copy_diag(data, data_tmp, dim1);
      }
      else
      {
         kernel_mult(data, data_tmp, dim1*dim2);
      }

      if(device != cu_M.device)
      {
         faust_cudaFree(data_tmp);
         data_tmp = NULL;
      }
   }

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_hadamard_product.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::Display()const
{
#ifdef __COMPILE_TIMERS__
t_display.start();
#endif
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M,*this);
   M.Display();
#ifdef __COMPILE_TIMERS__
t_display.stop();
#endif
}

template <typename FPP>
void Faust::MatDense<FPP,Gpu>::print_file(const char* filename)const
{
#ifdef __COMPILE_TIMERS__
t_print_file.start();
#endif
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M,*this);
   M.print_file(filename);
#ifdef __COMPILE_TIMERS__
t_print_file.stop();
#endif
}



#ifdef __COMPILE_TIMERS__
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_constructor_from_device;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_constructor_from_host;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_constructor_from_size;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_create;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_clear;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_copy_from_host;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_copy_from_device;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_copy_to_host;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_copy_to_device;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_init;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_move_to_device;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_set_zeros;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_set_eyes;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_transpose;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_init_from_transpose;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_max;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_min;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_abs;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_norm;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_trace;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_spectral_norm2;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_scalar_multiply;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_operator_equal_from_device;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_operator_equal_from_host;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_operator_plus_equal;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_operator_less_equal;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_add_cuspmat;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_sub;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_hadamard_product;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_display;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_print_file;

template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_gemm;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_gemv;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_add_ext;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_power_iteration;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_power_iteration_operator_equal;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_power_iteration_normalize;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_power_iteration_gemv;
template <typename FPP> Faust::Timer_gpu Faust::MatDense<FPP,Gpu>::t_power_iteration_dot;


template <typename FPP>
void Faust::MatDense<FPP,Gpu>::print_timers()const
{
   cout << "timers in Faust::MatDense :" << endl;
   cout << "t_constructor_from_device    = " << t_constructor_from_device.get_time()          << " s for "<< t_constructor_from_device.get_nb_call()          << " calls" << endl;
   cout << "t_constructor_from_host      = " << t_constructor_from_host.get_time()          << " s for "<< t_constructor_from_host.get_nb_call()          << " calls" << endl;
   cout << "t_constructor_from_size      = " << t_constructor_from_size.get_time()          << " s for "<< t_constructor_from_size.get_nb_call()          << " calls" << endl;
   cout << "t_create                     = " << t_create.get_time()          << " s for "<< t_create.get_nb_call()          << " calls" << endl;
   cout << "t_clear                      = " << t_clear.get_time()          << " s for "<< t_clear.get_nb_call()          << " calls" << endl;
   cout << "t_copy_from_host             = " << t_copy_from_host.get_time()          << " s for "<< t_copy_from_host.get_nb_call()          << " calls" << endl;
   cout << "t_copy_from_device           = " << t_copy_from_device.get_time()          << " s for "<< t_copy_from_device.get_nb_call()          << " calls" << endl;
   cout << "t_copy_to_host               = " << t_copy_to_host.get_time()          << " s for "<< t_copy_to_host.get_nb_call()          << " calls" << endl;
   cout << "t_copy_to_device             = " << t_copy_to_device.get_time()          << " s for "<< t_copy_to_device.get_nb_call()          << " calls" << endl;
   cout << "t_init                       = " << t_init.get_time()          << " s for "<< t_init.get_nb_call()          << " t_init" << endl;
   cout << "t_move_to_device             = " << t_move_to_device.get_time()          << " s for "<< t_move_to_device.get_nb_call()          << " calls" << endl;
   cout << "t_set_zeros                  = " << t_set_zeros.get_time()          << " s for "<< t_set_zeros.get_nb_call()          << " calls" << endl;
   cout << "t_set_eyes                   = " << t_set_eyes.get_time()          << " s for "<< t_set_eyes.get_nb_call()          << " calls" << endl;
   cout << "t_transpose                  = " << t_transpose.get_time()          << " s for "<< t_transpose.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_transpose        = " << t_init_from_transpose.get_time()          << " s for "<< t_init_from_transpose.get_nb_call()          << " calls" << endl;
   cout << "t_max                        = " << t_max.get_time()          << " s for "<< t_max.get_nb_call()          << " calls" << endl;
   cout << "t_min                        = " << t_min.get_time()          << " s for "<< t_min.get_nb_call()          << " calls" << endl;
   cout << "t_abs                        = " << t_abs.get_time()          << " s for "<< t_abs.get_nb_call()          << " calls" << endl;
   cout << "t_norm                       = " << t_norm.get_time()          << " s for "<< t_norm.get_nb_call()          << " calls" << endl;
   cout << "t_trace                      = " << t_trace.get_time()          << " s for "<< t_trace.get_nb_call()          << " calls" << endl;
   cout << "t_spectral_norm2             = " << t_spectral_norm2.get_time()          << " s for "<< t_spectral_norm2.get_nb_call()          << " calls" << endl;
   cout << "t_scalar_multiply            = " << t_scalar_multiply.get_time()          << " s for "<< t_scalar_multiply.get_nb_call()          << " calls" << endl;
   cout << "t_operator_equal_from_device = " << t_operator_equal_from_device.get_time()          << " s for "<< t_operator_equal_from_device.get_nb_call()          << " calls" << endl;
   cout << "t_operator_equal_from_host   = " << t_operator_equal_from_host.get_time()          << " s for "<< t_operator_equal_from_host.get_nb_call()          << " calls" << endl;
   cout << "t_operator_plus_equal        = " << t_operator_plus_equal.get_time()          << " s for "<< t_operator_plus_equal.get_nb_call()          << " calls" << endl;
   cout << "t_operator_less_equal        = " << t_operator_less_equal.get_time()          << " s for "<< t_operator_less_equal.get_nb_call()          << " calls" << endl;
   cout << "t_add_cuspmat                = " << t_add_cuspmat.get_time()          << " s for "<< t_add_cuspmat.get_nb_call()          << " calls" << endl;
   cout << "t_sub                        = " << t_sub.get_time()          << " s for "<< t_sub.get_nb_call()          << " calls" << endl;
   cout << "t_hadamard_product           = " << t_hadamard_product.get_time()          << " s for "<< t_hadamard_product.get_nb_call()          << " calls" << endl;
   cout << "t_display                    = " << t_display.get_time()          << " s for "<< t_display.get_nb_call()          << " calls" << endl;
   cout << "t_print_file                 = " << t_print_file.get_time()          << " s for "<< t_print_file.get_nb_call()          << " calls" << endl;
   cout << "t_gemm                       = " << t_gemm.get_time()          << " s for "<< t_gemm.get_nb_call()          << " calls" << endl;
   cout << "t_gemv                       = " << t_gemv.get_time()          << " s for "<< t_gemv.get_nb_call()          << " calls" << endl;
   cout << "t_add_ext                    = " << t_add_ext.get_time()          << " s for "<< t_add_ext.get_nb_call()          << " calls" << endl;
   cout << "t_power_iteration            = " << t_power_iteration.get_time()          << " s for "<< t_power_iteration.get_nb_call()          << " calls" << endl ;
   cout << "       t_power_iteration_operator_equal = " << t_power_iteration_operator_equal.get_time()          << " s for "<< t_power_iteration_operator_equal.get_nb_call()          << " calls" << endl;
   cout << "       t_power_iteration_normalize      = " << t_power_iteration_normalize.get_time()          << " s for "<< t_power_iteration_normalize.get_nb_call()          << " calls" << endl;
   cout << "       t_power_iteration_gemv           = " << t_power_iteration_gemv.get_time()          << " s for "<< t_power_iteration_gemv.get_nb_call()          << " calls" << endl ;
   cout << "       t_power_iteration_dot            = " << t_power_iteration_dot.get_time()          << " s for "<< t_power_iteration_dot.get_nb_call()          << " calls" << endl << endl <<endl;
}
#endif

#endif
