#ifndef __FAUST_CU_SPMAT_HPP__
#define __FAUST_CU_SPMAT_HPP__

#include "faust_MatDense.h"
#include "faust_MatDense_gpu.h"
#include "faust_MatSparse.h"
//#include <iostream>
#include <fstream>
#include <iomanip>
#include "faust_exception.h"
#include "kernels.h"
#include "faust_reduce_gpu.h"
#include "faust_gpu2cpu.h"
#include "faust_cuda.h"

using namespace std;

template <typename FPP>
const char * Faust::MatSparse<FPP,Gpu>::class_name="Faust::MatSparse<FPP,Gpu>::";

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse() :
	dim1(0),
	dim2(0),
	nnz(0),
	csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL){}

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse(const int* csrRowPtr_, const int* csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.start();
   else
      t_constructor_from_host.start();
#endif

    if(nnz_ == 0)
    {
        _create(nnz_, nbRow, nbCol, dstDevice);
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.stop();
   else
      t_constructor_from_host.stop();
#endif
        return;
    }
    if(csrRowPtr_==NULL || csrColInd_==NULL || csrValues_==NULL )
            handleError(class_name, "Faust::MatSparse<FPP,Gpu>::MatSparse(const faust_unsigned_int nnz_, const int* csrRowPtr_, const int* csrColInd_, const FPP* csrValues_, , const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE)");

    if(dataFromGPU)
    {
        copyFromDevice(csrRowPtr_, csrColInd_ ,csrValues_ , nnz, nbRow, nbCol, dstDevice, srcDevice, stream);
    }
   else
        copyFromHost(csrRowPtr_, csrColInd_, csrValues_, nnz, nbRow, nbCol, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
   if(dataFromGPU)
      t_constructor_from_device.stop();
   else
      t_constructor_from_host.stop();
#endif
}

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse(const Faust::MatSparse<FPP,Gpu>& cu_S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.start();
#endif
    _create(0,0,0);
    init(cu_S, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.stop();
#endif
}

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse(const Faust::MatSparse<FPP,Cpu>& S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.start();
#endif
    _create(0,0,0);
    init(S, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.stop();
#endif
}

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.start();
#endif
    _create(0,0,0);
    init(cu_A, spblasHandle, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_device.stop();
#endif
}

template <typename FPP>
Faust::MatSparse<FPP,Gpu>::MatSparse(const Faust::MatDense<FPP,Cpu>& A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.start();
#endif
    _create(0,0,0);
    init(A, spblasHandle, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_constructor_from_host.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::_create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const int device_)
{
#ifdef __COMPILE_TIMERS__
t_create.start();
#endif
    if (nnz_<0 || dim1_<0 || dim2_<0 || nnz_>dim1_*dim2_ )
       handleError(class_name, "_create : incorrect dimensions");



    if (dim1_ * dim2_ == 0)
    {
       if (nnz_ !=0 )
          handleError(class_name, "_create : to creat an empty matrix, nnz should be equal to 0");
       csrRowPtr = NULL;
       csrColInd = NULL;
       csrValues = NULL;
       descr = NULL;
    }
    else
    {
      int currentGPU;
       faust_cudaGetDevice(&currentGPU);
       faust_cudaSetDevice(device_);

       if (csrRowPtr == NULL)
          faust_cudaMalloc((void**)&csrRowPtr, (dim1_+1)*sizeof(int));
       else
          handleError(class_name, "_create : csrRowPtr has already been allocated on GPU");
       if(nnz_ > 0)
       {
          if (csrColInd == NULL)
             faust_cudaMalloc((void**)&csrColInd, (nnz_)*sizeof(int));
          else
             handleError(class_name, "_create : csrColInd has already been allocated on GPU");

          if (csrValues == NULL)
             faust_cudaMalloc((void**)&csrValues, (nnz_)*sizeof(FPP));
          else
             handleError(class_name, "_create : csrValues has already been allocated on GPU");
       }
       descr = &descr_content;
       faust_cusparseCreateMatDescr(descr);

       faust_cudaSetDevice(currentGPU);
    }

    nnz = nnz_;
    dim1 = dim1_;
    dim2 = dim2_;
    device = device_;
#ifdef __COMPILE_TIMERS__
t_create.stop();
#endif
}


template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::_clear()
{
#ifdef __COMPILE_TIMERS__
t_clear.start();
#endif
    if(dim1*dim2 > 0)
    {
        int currentGPU;
        faust_cudaGetDevice(&currentGPU);
        faust_cudaSetDevice(device);

        if (csrRowPtr!=NULL)
           faust_cudaFree(csrRowPtr);
        else
           handleError(class_name, "_clear : csrRowPtr has not been allocated");

        if(nnz > 0)
        {
           if (csrColInd!=NULL)
              faust_cudaFree(csrColInd);
           else
              handleError(class_name, "_clear : csrColInd has not been allocated");

           if (csrValues!=NULL)
              faust_cudaFree(csrValues);
           else
              handleError(class_name, "_clear : csrValues has not been allocated");
       }
       faust_cusparseDestroyMatDescr(*descr);
       faust_cudaSetDevice(currentGPU);
    }
    else if(descr != NULL)
       handleError(class_name, "_clear : descr pointer for empty sparse matrix should be NULL");

    nnz = 0;
    dim1 = 0;
    dim2 = 0;
    csrRowPtr = NULL;
    csrColInd = NULL;
    csrValues = NULL;
    descr = NULL;

#ifdef __COMPILE_TIMERS__
t_clear.stop();
#endif

}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_)
{
   if (nnz_!=nnz || dim1!=nbRow || device!=device_ )
   {
      _clear();
      _create(nnz_, nbRow, nbCol, device_);
   }
   else
   {
      nnz = nnz_;
      dim1 = nbRow;
      dim2 = nbCol;
      device = device_;
   }

}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::copyFromHost(const int* csrRowPtr_, const int*  csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_copy_from_host.start();
#endif
    resize(nnz_, nbRow, nbCol, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(csrRowPtr, csrRowPtr_, (nbRow+1)*sizeof(int), cudaMemcpyHostToDevice, stream);
    faust_cudaMemcpyAsync(csrColInd, csrColInd_, nnz_*sizeof(int), cudaMemcpyHostToDevice, stream);
    faust_cudaMemcpyAsync(csrValues, csrValues_, nnz_*sizeof(FPP), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_copy_from_host.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::copyFromDevice(const int* csrRowPtr_, const int* csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_,  const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
#ifdef __COMPILE_TIMERS__
t_copy_from_device.start();
#endif
    resize(nnz_, nbRow, nbCol, dstDevice);
    faust_cudaMemcpyPeerAsync(csrRowPtr, dstDevice, csrRowPtr_, srcDevice, (nbRow+1)*sizeof(int), stream);
    faust_cudaMemcpyPeerAsync(csrColInd, dstDevice, csrColInd_, srcDevice, nnz_*sizeof(int), stream);
    faust_cudaMemcpyPeerAsync(csrValues, dstDevice, csrValues_, srcDevice, nnz_*sizeof(FPP), stream);
#ifdef __COMPILE_TIMERS__
t_copy_from_device.stop();
#endif
}


template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::copyToHost(int* csrRowPtr_, int* csrColInd_, FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream/*=0*/)const
{
#ifdef __COMPILE_TIMERS__
t_copy_to_host.start();
#endif
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(csrRowPtr_, csrRowPtr, (nbRow+1)*sizeof(int), cudaMemcpyDeviceToHost, stream);
    faust_cudaMemcpyAsync(csrColInd_, csrColInd, nnz_*sizeof(int), cudaMemcpyDeviceToHost, stream);
    faust_cudaMemcpyAsync(csrValues_, csrValues, nnz_*sizeof(FPP), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_copy_to_host.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::copyToDevice(int* csrRowPtr_, int* csrColInd_, FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
#ifdef __COMPILE_TIMERS__
t_copy_to_device.start();
#endif
       faust_cudaMemcpyPeerAsync(csrRowPtr_, dstDevice, csrRowPtr, device, (nbRow+1)*sizeof(int), stream);
       faust_cudaMemcpyPeerAsync(csrColInd_, dstDevice, csrColInd, device, nnz_*sizeof(int), stream);
       faust_cudaMemcpyPeerAsync(csrValues_, dstDevice, csrValues, device, nnz_*sizeof(FPP), stream);
#ifdef __COMPILE_TIMERS__
t_copy_to_device.stop();
#endif
}



template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
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
   if (nnz == 0)
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

   int* rowptr_tmp;
   int* colind_tmp;
   FPP* values_tmp;
   faust_cudaMalloc((void**)&rowptr_tmp, (dim1+1)*sizeof(int));
   faust_cudaMalloc((void**)&colind_tmp, nnz*sizeof(int));
   faust_cudaMalloc((void**)&values_tmp, nnz*sizeof(FPP));
   faust_cudaMemcpyPeerAsync(rowptr_tmp, dstDevice, csrRowPtr, device, (dim1+1)*sizeof(int), stream);
   faust_cudaMemcpyPeerAsync(colind_tmp, dstDevice, csrColInd, device, nnz*sizeof(int), stream);
   faust_cudaMemcpyPeerAsync(values_tmp, dstDevice, csrValues, device, nnz*sizeof(FPP), stream);
   faust_cudaFree(csrRowPtr);
   faust_cudaFree(csrColInd);
   faust_cudaFree(csrValues);
   csrRowPtr = rowptr_tmp;
   csrColInd = colind_tmp;
   csrValues = values_tmp;
   device = dstDevice;

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_move_to_device.stop();
#endif
}


template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::init(const Faust::MatSparse<FPP,Gpu>& cu_S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
#ifdef __COMPILE_TIMERS__
t_init_from_cuspmat.start();
#endif

    if(cu_S.nnz == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.nnz, cu_S.dim1, cu_S.dim2);
       device = cu_S.device;
#ifdef __COMPILE_TIMERS__
t_init_from_cuspmat.stop();
#endif
       return;
    }
    else
       copyFromDevice(cu_S.csrRowPtr, cu_S.csrColInd, cu_S.csrValues, cu_S.nnz, cu_S.getNbRow(), cu_S.getNbCol(), cu_S.device, cu_S.device);
#ifdef __COMPILE_TIMERS__
t_init_from_cuspmat.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::init(const Faust::MatSparse<FPP,Cpu>& S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
#ifdef __COMPILE_TIMERS__
t_init_from_spmat.start();
#endif
    if(S.getNonZeros() == 0)
        resize(S.getNonZeros(), S.getNbRow(), S.getNbCol(), dstDevice);
    else
        copyFromHost(S.getRowPtr(), S.getColInd(), S.getValuePtr(), S.getNonZeros(), S.getNbRow(), S.getNbCol(), dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_init_from_spmat.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::init(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
#ifdef __COMPILE_TIMERS__
t_init_from_cumat.start();
#endif
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    const Faust::MatDense<FPP,Gpu>* cu_A_ptr;

    if(device != dstDevice)
       cu_A_ptr = new Faust::MatDense<FPP,Gpu>(cu_A, dstDevice);
    else
       cu_A_ptr = &cu_A;

    int dim1_= (int)cu_A_ptr->getNbRow();
    int dim2_= (int)cu_A_ptr->getNbCol();
    int nnz_=0;



    if(cu_A_ptr->getData())
    {
       int* nnzPerRow;
       faust_cudaMalloc((void**)&nnzPerRow, dim1_*sizeof(int));

       cusparseMatDescr_t descr_tmp;
       faust_cusparseCreateMatDescr(&descr_tmp);

       faust_cu_nnz(spblasHandle.GetHandle(), CUSPARSE_DIRECTION_ROW,
             dim1_, dim2_, descr_tmp,
             cu_A_ptr->getData(),
             dim1_, nnzPerRow, &nnz_);

       resize(nnz_, dim1_, dim2_, dstDevice);

       faust_cu_dense2csr(spblasHandle.GetHandle(),
             dim1_, dim2_, *descr,
             cu_A_ptr->getData(),
             dim1_, nnzPerRow,
             csrValues, csrRowPtr, csrColInd);
       faust_cudaFree(nnzPerRow);
    }
    else if(cu_A_ptr->estIdentite())
    {
       nnz_ = dim1_;
       resize(nnz_, dim1_, dim1_, dstDevice);
       int* rowptr = new int[nnz_+1];
       for(int i=0 ; i<nnz_+1 ; i++)
          rowptr[i] =  i;
       int* colind = new int[nnz_];
       for(int i=0 ; i<nnz_ ; i++)
          colind[i] = i;
       FPP* values = new FPP[nnz_];
       for(int i=0 ; i<nnz_ ; i++)
          values[i] = 1.0;

       copyFromHost(rowptr, colind, values, nnz_, dim1_, dim2_, dstDevice);
    }
    else if(cu_A_ptr->estNulle())
    {
       if(device != dstDevice)
          delete cu_A_ptr;
       handleError(class_name, "init(faust_cu_mat, ...) : zero sparse matrix is not managed. Please keep full matrix representation");
    }
    else
    {
       if(device != dstDevice)
          delete cu_A_ptr;
       handleError(class_name, "init(faust_cu_mat, ...) : impossible to create sparse matrix from uninitialized full matrix");
    }

    if(device != dstDevice)
       delete cu_A_ptr;
    cu_A_ptr = NULL;

    faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_init_from_cumat.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::init(const Faust::MatDense<FPP,Cpu>& M, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
#ifdef __COMPILE_TIMERS__
t_init_from_mat.start();
#endif
    const Faust::MatDense<FPP,Gpu> cu_M(M.getNbRow(), M.getNbCol(), dstDevice);
    init(cu_M, spblasHandle, dstDevice, stream);
#ifdef __COMPILE_TIMERS__
t_init_from_mat.stop();
#endif
}



template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::operator*=(const FPP alpha)
{
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.start();
#endif
   if (alpha == 1.0 || nnz == 0)
   {
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
      return;
   }
   if (alpha == 0.0)
      resize(0, dim1, dim2);
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_mult_const(csrValues, alpha, nnz);

      faust_cudaSetDevice(currentGPU);
   }
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
}

template <typename FPP>
bool Faust::MatSparse<FPP,Gpu>::operator==(const Faust::MatSparse<FPP,Gpu>& cu_S)const
{
   if(csrRowPtr!=NULL && csrRowPtr==cu_S.csrRowPtr && device==cu_S.device && nnz==cu_S.nnz && csrColInd==cu_S.csrColInd && csrValues==cu_S.csrValues && dim2==cu_S.dim2)
   {
      if(dim1!=cu_S.dim1 || descr!=cu_S.descr)
         handleError(class_name,"operator== : same data and device but different dimensions");
      return true;
   }
   else
      return false;
}
template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::operator+=(const FPP alpha)
{
#ifdef __COMPILE_TIMERS__
t_operator_plus_equal_real.start();
#endif
   if (alpha == 0.0 || nnz == 0)
   {
#ifdef __COMPILE_TIMERS__
t_operator_plus_equal_real.stop();
#endif
      return;
   }
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_add_const(csrValues, alpha, nnz);

      faust_cudaSetDevice(currentGPU);
   }
#ifdef __COMPILE_TIMERS__
t_operator_plus_equal_real.stop();
#endif
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::init_from_transpose(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.start();
#endif
   resize(cu_S.nnz, cu_S.dim2, cu_S.dim1);

   if(cu_S.dim1*cu_S.dim2==0 || cu_S.nnz==0)
   {
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif
      return;
   }


   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);


  if(device != cu_S.device)
  {
      int *rowPtrOld, *colIndOld;
      faust_cudaMalloc((void**)&rowPtrOld, (cu_S.dim1+1)*sizeof(int));
      faust_cudaMalloc((void**)&colIndOld, cu_S.nnz*sizeof(int));
      faust_cudaMemcpy(rowPtrOld, cu_S.getRowPtr(), (cu_S.dim1+1)*sizeof(int), cudaMemcpyDeviceToDevice);
      faust_cudaMemcpy(colIndOld, cu_S.getColInd(), cu_S.nnz*sizeof(int), cudaMemcpyDeviceToDevice);

      faust_cu_csr2csc(spblasHandle.GetHandle(),
         dim1, dim2, nnz,
         (FPP*)NULL, rowPtrOld, colIndOld,
         (FPP*)NULL, csrColInd, csrRowPtr,
         CUSPARSE_ACTION_SYMBOLIC,
         CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaFree(rowPtrOld);
      faust_cudaFree(colIndOld);

   }
   else
   {
      faust_cu_csr2csc(spblasHandle.GetHandle(),
         dim1, dim2, nnz,
         (FPP*)NULL, cu_S.getRowPtr(), cu_S.getColInd(),
         (FPP*)NULL, csrColInd, csrRowPtr,
         CUSPARSE_ACTION_SYMBOLIC,
         CUSPARSE_INDEX_BASE_ZERO);

   }


   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_init_from_transpose.stop();
#endif
}
template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::transpose(Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif
   if(dim1*dim2==0 || nnz==0)
   {
      resize(nnz,dim2,dim1);
      return;
   }

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(device);

   int *rowPtrOld, *colIndOld;
   faust_cudaMalloc((void**)&rowPtrOld, (dim1+1)*sizeof(int));
   faust_cudaMalloc((void**)&colIndOld, nnz*sizeof(int));
   faust_cudaMemcpy(rowPtrOld, getRowPtr(), (dim1+1)*sizeof(int), cudaMemcpyDeviceToDevice);
   faust_cudaMemcpy(colIndOld, getColInd(), nnz*sizeof(int), cudaMemcpyDeviceToDevice);

   resize(nnz, dim2, dim1);

   faust_cu_csr2csc(spblasHandle.GetHandle(),
         dim1, dim2, nnz,
         (FPP*)NULL, rowPtrOld, colIndOld,
         (FPP*)NULL, csrColInd, csrRowPtr,
         CUSPARSE_ACTION_SYMBOLIC,
         CUSPARSE_INDEX_BASE_ZERO);

   faust_cudaFree(rowPtrOld);
   faust_cudaFree(colIndOld);

   faust_cudaSetDevice(currentGPU);
#ifdef __COMPILE_TIMERS__
t_transpose.stop();
#endif
}

template <typename FPP>
FPP Faust::MatSparse<FPP,Gpu>::norm() const
{
#ifdef __COMPILE_TIMERS__
t_norm.start();
#endif
   if (nnz>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      FPP val =  faust_cu_norm(csrValues, nnz);
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
}

template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::Display()const
{
#ifdef __COMPILE_TIMERS__
t_display.start();
#endif
   Faust::MatSparse<FPP,Cpu> S;
   faust_gpu2cpu(S,*this);
   S.Display();
#ifdef __COMPILE_TIMERS__
t_display.stop();
#endif
}

template<typename FPP>
void Faust::MatSparse<FPP,Gpu>::print_file(const char* filename)const
{print_file(filename,std::fstream::out);}

template<typename FPP>
void Faust::MatSparse<FPP,Gpu>::print_file(const char* filename,std::ios_base::openmode mode)const
{
#ifdef __COMPILE_TIMERS__
t_print_file.start();
#endif
   Faust::MatSparse<FPP,Cpu> tmp;
   faust_gpu2cpu(tmp, *this);
   tmp.print_file(filename, mode);
#ifdef __COMPILE_TIMERS__
t_print_file.stop();
#endif
}

#ifdef __COMPILE_TIMERS__

template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_constructor_from_device;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_constructor_from_host;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_create;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_clear;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_copy_from_host;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_copy_from_device;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_copy_to_host;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_copy_to_device;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_move_to_device;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_init_from_cuspmat;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_init_from_spmat;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_init_from_cumat;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_init_from_mat;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_scalar_multiply;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_operator_plus_equal_real;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_transpose;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_init_from_transpose;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_norm;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_display;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_print_file;

template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_csrmv;
template <typename FPP> Faust::Timer_gpu Faust::MatSparse<FPP,Gpu>::t_csrmm;


template <typename FPP>
void Faust::MatSparse<FPP,Gpu>::print_timers()const
{

   cout << "timers in faust_cu_spmat :" << endl;
   cout << "t_constructor_from_device   = " << t_constructor_from_device.get_time()          << " s for "<< t_constructor_from_device.get_nb_call()          << " calls" << endl;
   cout << "t_constructor_from_host     = " << t_constructor_from_host.get_time()          << " s for "<< t_constructor_from_host.get_nb_call()          << " calls" << endl;
   cout << "t_create                    = " << t_create.get_time()          << " s for "<< t_create.get_nb_call()          << " calls" << endl;
   cout << "t_clear                     = " << t_clear.get_time()          << " s for "<< t_clear.get_nb_call()          << " calls" << endl;
   cout << "t_copy_from_host            = " << t_copy_from_host.get_time()          << " s for "<< t_copy_from_host.get_nb_call()          << " calls" << endl;
   cout << "t_copy_from_device          = " << t_copy_from_device.get_time()          << " s for "<< t_copy_from_device.get_nb_call()          << " calls" << endl;
   cout << "t_copy_to_host              = " << t_copy_to_host.get_time()          << " s for "<< t_copy_to_host.get_nb_call()          << " calls" << endl;
   cout << "t_copy_to_device            = " << t_copy_to_device.get_time()          << " s for "<< t_copy_to_device.get_nb_call()          << " calls" << endl;
   cout << "t_move_to_device            = " << t_move_to_device.get_time()          << " s for "<< t_move_to_device.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_cuspmat         = " << t_init_from_cuspmat.get_time()          << " s for "<< t_init_from_cuspmat.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_spmat           = " << t_init_from_spmat.get_time()          << " s for "<< t_init_from_spmat.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_cumat           = " << t_init_from_cumat.get_time()          << " s for "<< t_init_from_cumat.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_mat             = " << t_init_from_mat.get_time()          << " s for "<< t_init_from_mat.get_nb_call()          << " calls" << endl;
   cout << "t_scalar_multiply           = " << t_scalar_multiply.get_time()          << " s for "<< t_scalar_multiply.get_nb_call()          << " calls" << endl;
   cout << "t_operator_plus_equal_real  = " << t_operator_plus_equal_real.get_time()          << " s for "<< t_operator_plus_equal_real.get_nb_call()          << " calls" << endl;
   cout << "t_transpose                 = " << t_transpose.get_time()          << " s for "<< t_transpose.get_nb_call()          << " calls" << endl;
   cout << "t_init_from_transpose       = " << t_init_from_transpose.get_time()          << " s for "<< t_init_from_transpose.get_nb_call()          << " calls" << endl;
   cout << "t_norm                      = " << t_norm.get_time()          << " s for "<< t_norm.get_nb_call()          << " calls" << endl;
   cout << "t_display                   = " << t_display.get_time()          << " s for "<< t_display.get_nb_call()          << " calls" << endl;
   cout << "t_print_file                = " << t_print_file.get_time()          << " s for "<< t_print_file.get_nb_call()          << " calls" << endl;
   cout << "t_csrmv                     = " << t_csrmv.get_time()          << " s for "<< t_csrmv.get_nb_call()          << " calls" << endl;
   cout << "t_csrmm                     = " << t_csrmm.get_time()          << " s for "<< t_csrmm.get_nb_call()          << " calls" << endl << endl << endl;

}

#endif

#endif
