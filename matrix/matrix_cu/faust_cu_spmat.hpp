#ifndef __FAUST_CU_SPMAT_HPP__
#define __FAUST_CU_SPMAT_HPP__

#include "faust_mat.h"
#include "faust_cu_mat.h"
#include "faust_spmat.h"
//#include <iostream>
#include <fstream>
#include <iomanip>
#include "faust_exception.h"
#include "kernels.h"
#include "faust_cu_reduce.h"
#include "faust_cu2faust.h"

using namespace std;

template <typename faust_real>
const char * faust_cu_spmat<faust_real>::class_name="faust_cu_spmat<faust_real>::";

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat() :
	dim1(0),
	dim2(0),
	nnz(0),
	csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL){}

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat(const int* csrRowPtr_, const int* csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
    if(nnz_ == 0)
    {
        _create(nnz_, nbRow, nbCol, dstDevice);
        return;
    }
    if(csrRowPtr_==NULL || csrColInd_==NULL || csrValues_==NULL )
            handleError(class_name, "faust_cu_spmat<faust_real>::faust_cu_spmat(const faust_unsigned_int nnz_, const int* csrRowPtr_, const int* csrColInd_, const faust_real* csrValues_, , const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL),device(FAUST_DEFAULT_CUDA_DEVICE)");

    if(dataFromGPU)
    {
        copyFromDevice(csrRowPtr_, csrColInd_ ,csrValues_ , nnz, nbRow, nbCol, dstDevice, srcDevice, stream);
    }
   else
        copyFromHost(csrRowPtr_, csrColInd_, csrValues_, nnz, nbRow, nbCol, dstDevice, stream);
}

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat(const faust_cu_spmat<faust_real>& cu_S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
    _create(0,0,0);
    init(cu_S, dstDevice, stream);
}

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat(const faust_spmat<faust_real>& S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
    _create(0,0,0);
    init(S, dstDevice, stream);
}

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat(const faust_cu_mat<faust_real>& cu_A, cusparseHandle_t cusparseHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
    _create(0,0,0);
    init(cu_A, cusparseHandle, dstDevice, stream);
}

template <typename faust_real>
faust_cu_spmat<faust_real>::faust_cu_spmat(const faust_mat<faust_real>& A, cusparseHandle_t cusparseHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), nnz(0), csrRowPtr(NULL), csrColInd(NULL), csrValues(NULL), device(FAUST_DEFAULT_CUDA_DEVICE), descr(NULL)
{
    _create(0,0,0);
    init(A, cusparseHandle, dstDevice, stream);
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::_create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const int device_)
{
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
             faust_cudaMalloc((void**)&csrValues, (nnz_)*sizeof(faust_real));
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
}


template <typename faust_real>
void faust_cu_spmat<faust_real>::_clear()
{
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


}

template <typename faust_real>
void faust_cu_spmat<faust_real>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_)
{
#ifdef __COMPILE_TIMERS__
t_resize.start();
#endif
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

#ifdef __COMPILE_TIMERS__
t_resize.stop();
#endif
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::copyFromHost(const int* csrRowPtr_, const int*  csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    resize(nnz_, nbRow, nbCol, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(csrRowPtr, csrRowPtr_, (nbRow+1)*sizeof(int), cudaMemcpyHostToDevice, stream);
    faust_cudaMemcpyAsync(csrColInd, csrColInd_, nnz_*sizeof(int), cudaMemcpyHostToDevice, stream);
    faust_cudaMemcpyAsync(csrValues, csrValues_, nnz_*sizeof(faust_real), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::copyFromDevice(const int* csrRowPtr_, const int* csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_,  const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    resize(nnz_, nbRow, nbCol, dstDevice);
    faust_cudaMemcpyPeerAsync(csrRowPtr, dstDevice, csrRowPtr_, srcDevice, (nbRow+1)*sizeof(int), stream);
    faust_cudaMemcpyPeerAsync(csrColInd, dstDevice, csrColInd_, srcDevice, nnz_*sizeof(int), stream);
    faust_cudaMemcpyPeerAsync(csrValues, dstDevice, csrValues_, srcDevice, nnz_*sizeof(faust_real), stream);
}


template <typename faust_real>
void faust_cu_spmat<faust_real>::copyToHost(int* csrRowPtr_, int* csrColInd_, faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream/*=0*/)const
{
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(csrRowPtr_, csrRowPtr, (nbRow+1)*sizeof(int), cudaMemcpyDeviceToHost, stream);
    faust_cudaMemcpyAsync(csrColInd_, csrColInd, nnz_*sizeof(int), cudaMemcpyDeviceToHost, stream);
    faust_cudaMemcpyAsync(csrValues_, csrValues, nnz_*sizeof(faust_real), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::copyToDevice(int* csrRowPtr_, int* csrColInd_, faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
       faust_cudaMemcpyPeerAsync(csrRowPtr_, dstDevice, csrRowPtr, device, (nbRow+1)*sizeof(int), stream);
       faust_cudaMemcpyPeerAsync(csrColInd_, dstDevice, csrColInd, device, nnz_*sizeof(int), stream);
       faust_cudaMemcpyPeerAsync(csrValues_, dstDevice, csrValues, device, nnz_*sizeof(faust_real), stream);
}



template <typename faust_real>
void faust_cu_spmat<faust_real>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
   if(device == dstDevice)
      return;
   if (nnz == 0)
   {
      device = dstDevice;
      return;
   }

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(dstDevice);

   int* rowptr_tmp;
   int* colind_tmp;
   faust_real* values_tmp;
   faust_cudaMalloc((void**)&rowptr_tmp, (dim1+1)*sizeof(int));
   faust_cudaMalloc((void**)&colind_tmp, nnz*sizeof(int));
   faust_cudaMalloc((void**)&values_tmp, nnz*sizeof(faust_real));
   faust_cudaMemcpyPeerAsync(rowptr_tmp, dstDevice, csrRowPtr, device, (dim1+1)*sizeof(int), stream);
   faust_cudaMemcpyPeerAsync(colind_tmp, dstDevice, csrColInd, device, nnz*sizeof(int), stream);
   faust_cudaMemcpyPeerAsync(values_tmp, dstDevice, csrValues, device, nnz*sizeof(faust_real), stream);
   faust_cudaFree(csrRowPtr);
   faust_cudaFree(csrColInd);
   faust_cudaFree(csrValues);
   csrRowPtr = rowptr_tmp;
   csrColInd = colind_tmp;
   csrValues = values_tmp;
   device = dstDevice;

   faust_cudaSetDevice(currentGPU);
}


template <typename faust_real>
void faust_cu_spmat<faust_real>::init(const faust_cu_spmat<faust_real>& cu_S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{

    if(cu_S.nnz == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.nnz, cu_S.dim1, cu_S.dim2);
       device = cu_S.device;
       return;
    }
    else
       copyFromDevice(cu_S.csrRowPtr, cu_S.csrColInd, cu_S.csrValues, cu_S.nnz, cu_S.getNbRow(), cu_S.getNbCol(), cu_S.device, cu_S.device);
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::init(const faust_spmat<faust_real>& S, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
    if(S.getNonZeros() == 0)
        resize(S.getNonZeros(), S.getNbRow(), S.getNbCol(), dstDevice);
    else
        copyFromHost(S.getRowPtr(), S.getColInd(), S.getValuePtr(), S.getNonZeros(), S.getNbRow(), S.getNbCol(), dstDevice, stream);
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::init(const faust_cu_mat<faust_real>& cu_A, cusparseHandle_t cusparseHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    const faust_cu_mat<faust_real>* cu_A_ptr;

    if(device != dstDevice)
       cu_A_ptr = new faust_cu_mat<faust_real>(cu_A, dstDevice);
    else
       cu_A_ptr = &cu_A;

    const int dim1_= (int)cu_A_ptr->getNbRow();
    const int dim2_= (int)cu_A_ptr->getNbCol();
    int nnz_;

    if(cu_A_ptr->getData())
    {
       int* nnzPerRow;
       faust_cudaMalloc((void**)&nnzPerRow, dim1_*sizeof(int));
       faust_cu_nnz(cusparseHandle, CUSPARSE_DIRECTION_ROW,
             dim1_, dim2_, *descr,
             cu_A_ptr->getData(),
             dim1_, nnzPerRow, &nnz_);

       resize(nnz_, dim1_, dim2_, dstDevice);

       faust_cu_dense2csr(cusparseHandle,
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
       faust_real* values = new faust_real[nnz_];
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
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::init(const faust_mat<faust_real>& M, cusparseHandle_t cusparseHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ )
{
    const faust_cu_mat<faust_real> cu_M(M.getNbRow(), M.getNbCol(), dstDevice);
    init(cu_M, cusparseHandle, dstDevice, stream);
}



template <typename faust_real>
void faust_cu_spmat<faust_real>::operator*=(const faust_real alpha)
{
   if (alpha == 1.0 || nnz == 0)
      return;
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
}

template <typename faust_real>
bool faust_cu_spmat<faust_real>::operator==(const faust_cu_spmat<faust_real>& cu_S)const
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
template <typename faust_real>
void faust_cu_spmat<faust_real>::operator+=(const faust_real alpha)
{
   if (alpha == 0.0 || nnz == 0)
      return;
   else
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      kernel_add_const(csrValues, alpha, nnz);

      faust_cudaSetDevice(currentGPU);
   }
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::init_from_transpose(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle)
{
   resize(cu_S.nnz, cu_S.dim2, cu_S.dim1);

   if(cu_S.dim1*cu_S.dim2==0 || cu_S.nnz==0)
      return;
   

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

      faust_cu_csr2csc(cusparseHandle,
         dim1, dim2, nnz,
         (faust_real*)NULL, rowPtrOld, colIndOld, 
         (faust_real*)NULL, csrColInd, csrRowPtr, 
         CUSPARSE_ACTION_SYMBOLIC, 
         CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaFree(rowPtrOld);
      faust_cudaFree(colIndOld);

   }   
   else
   {
      faust_cu_csr2csc(cusparseHandle,
         dim1, dim2, nnz,
         (faust_real*)NULL, cu_S.getRowPtr(), cu_S.getColInd(), 
         (faust_real*)NULL, csrColInd, csrRowPtr, 
         CUSPARSE_ACTION_SYMBOLIC, 
         CUSPARSE_INDEX_BASE_ZERO);
      
   }


   faust_cudaSetDevice(currentGPU); 
}
template <typename faust_real>
void faust_cu_spmat<faust_real>::transpose(cusparseHandle_t cusparseHandle)
{
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

   faust_cu_csr2csc(cusparseHandle,
         dim1, dim2, nnz,
         (faust_real*)NULL, rowPtrOld, colIndOld, 
         (faust_real*)NULL, csrColInd, csrRowPtr, 
         CUSPARSE_ACTION_SYMBOLIC, 
         CUSPARSE_INDEX_BASE_ZERO);
   
   faust_cudaFree(rowPtrOld);
   faust_cudaFree(colIndOld);

   faust_cudaSetDevice(currentGPU); 
}

template <typename faust_real>
faust_real faust_cu_spmat<faust_real>::norm() const 
{
   if (nnz>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val =  faust_cu_norm(csrValues, nnz);
      faust_cudaSetDevice(currentGPU);
      return val;
   }
   else
   {
      // display warning as matrix is empty
      return (faust_real)0.0;
   }
}

template <typename faust_real>
void faust_cu_spmat<faust_real>::Display()const
{
   faust_spmat<faust_real> S;
   faust_cu2faust(S,*this);
   S.Display();
}

template<typename faust_real> 
void faust_cu_spmat<faust_real>::print_file(const char* filename)const
{print_file(filename,std::fstream::out);}

template<typename faust_real> 
void faust_cu_spmat<faust_real>::print_file(const char* filename,std::ios_base::openmode mode)const
{
   faust_spmat<faust_real> tmp;
   faust_cu2faust(tmp, *this);
   tmp.print_file(filename, mode);
}

#endif
