#ifndef __FAUST_CU_MAT_HPP__
#define __FAUST_CU_MAT_HPP__

//#include "faust_cu_mat.h"
#include "faust_cu_vec.h"
#include "faust_mat.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifdef __COMPILE_SPMAT__
  #include "faust_spmat.h"
  #include "faust_cu_spmat.h"
#endif
#include "faust_cu_reduce.h"
#include "kernels.h"
#include "faust_cuda.h"
#include "faust_cu2faust.h"

#include "LinAlgebra_cu.h"


using namespace std;

template <typename faust_real>
const char * faust_cu_mat<faust_real>::class_name = "faust_cu_mat<faust_real>::";

template <typename faust_real>
faust_cu_mat<faust_real>::faust_cu_mat(const faust_real *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    if(nbRow*nbCol == 0)
    {
        resize(nbRow, nbCol, dstDevice);
        return;
    }
    if(data_ == NULL)
            handleError(class_name, "faust_cu_mat(const faust_real*,const faust_unsigned_int, const faust_unsigned_int, bool, int=FAUST_DEFAULT_CUDA_DEVICE, int=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t=0) : data pointer argument is NULL");

    if(dataFromGPU)
        copyFromDevice(data_, nbRow, nbCol, dstDevice, srcDevice, stream);
    else
        copyFromHost(data_, nbRow, nbCol, dstDevice, stream);
}



template <typename faust_real>
faust_cu_mat<faust_real>::faust_cu_mat(const faust_cu_mat<faust_real>& cu_M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
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
}

template <typename faust_real>
faust_cu_mat<faust_real>::faust_cu_mat(const faust_mat<faust_real>& M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream /*=0*/ ) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL),device(FAUST_DEFAULT_CUDA_DEVICE)
{
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
}

#ifdef __COMPILE_SPMAT__
template <typename faust_real>
faust_cu_mat<faust_real>::faust_cu_mat(const faust_cu_spmat<faust_real>& cu_S,cusparseHandle_t cusparseHandle, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    if(cu_S.getNonZeros() == 0)
    {
        if(cu_S.getRowPtr()!=NULL || cu_S.getColInd()!=NULL || cu_S.getValues()!=NULL)
            handleError(class_name, "faust_cu_mat(const faust_cu_spmat&, ...) : rowPtr,  colInd or values pointer is not NULL");

        setZeros(cu_S.getNbRow(), cu_S.getNbCol(), dstDevice);
        return;
    }

    resize(cu_S.getNbRow(), cu_S.getNbCol(), dstDevice);

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    const faust_cu_spmat<faust_real>* cu_S_dst = &cu_S;
    if(dstDevice != cu_S.getDevice())
        cu_S_dst = new faust_cu_spmat<faust_real>(cu_S, dstDevice);

    faust_cu_csr2dense(cusparseHandle,
            dim1,dim2,
            cu_S_dst->getDescr(),
            cu_S_dst->getValues(), cu_S_dst->getRowPtr(), cu_S_dst->getColInd(),
            data, dim1);


    if(dstDevice != cu_S.getDevice())
        delete    cu_S_dst;
    cu_S_dst = NULL;

   faust_cudaSetDevice(currentGPU);
}
#endif

template <typename faust_real>
faust_cu_mat<faust_real>::faust_cu_mat(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/) : dim1(0), dim2(0), isIdentity(false), isZeros(false), data(NULL), device(FAUST_DEFAULT_CUDA_DEVICE)
{
    resize(nbRow, nbCol, dstDevice);
}



/// GETTEUR SETTEUR ///

template <typename faust_real>
void faust_cu_mat<faust_real>::_create(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_)
{
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
         faust_cudaMalloc((void**)&data, (nbRow*nbCol)*sizeof(faust_real));
      else
         handleError(class_name, "_create : data has already been allocated on GPU");
         
      
      faust_cudaSetDevice(currentGPU);
   }
   dim1 = nbRow;
   dim2 = nbCol;
   device = device_;
   isZeros = false;
   isIdentity = false;

}

template <typename faust_real>
void faust_cu_mat<faust_real>::_clear()
{
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
{
      handleError(class_name, "_clear : data of empty matrix, identity matrix or zeros matrix should be NULL");
}
   dim1 = 0;
   dim2 = 0;
   isIdentity = false;
   isZeros = false;
   data = NULL;
}


template <typename faust_real>
void faust_cu_mat<faust_real>::resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_)
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

template <typename faust_real>
void faust_cu_mat<faust_real>::copyFromHost(const faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    resize(nbRow, nbCol, dstDevice);
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(dstDevice);

    faust_cudaMemcpyAsync(data, data_, nbRow*nbCol*sizeof(faust_real), cudaMemcpyHostToDevice, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_mat<faust_real>::copyFromDevice(const faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, int srcDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
    resize(nbRow, nbCol, dstDevice);
    faust_cudaMemcpyPeerAsync(data, dstDevice, data_, srcDevice, nbRow*nbCol*sizeof(faust_real), stream);
}

template <typename faust_real>
void faust_cu_mat<faust_real>::copyToHost(faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream/*=0*/)const
{
    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpyAsync(data_, data, nbRow*nbCol*sizeof(faust_real), cudaMemcpyDeviceToHost, stream);
    faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_mat<faust_real>::copyToDevice(faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)const
{
       faust_cudaMemcpyPeerAsync(data_, dstDevice, data, device, nbRow*nbCol*sizeof(faust_real), stream);
}

template <typename faust_real>
void faust_cu_mat<faust_real>::init(const faust_cu_mat<faust_real>& cu_M, int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
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

}

template <typename faust_real>
void faust_cu_mat<faust_real>::moveToDevice(int dstDevice/*=FAUST_DEFAULT_CUDA_DEVICE*/, cudaStream_t stream/*=0*/)
{
   if(device == dstDevice)
      return;
   if (data == NULL)
   {
      device = dstDevice;
      return;
   }

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   faust_cudaSetDevice(dstDevice);

   faust_real* data_tmp;
   faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(faust_real));
   faust_cudaMemcpyPeerAsync(data_tmp, dstDevice, data, device, dim1*dim2*sizeof(faust_real), stream);
   faust_cudaFree(data);
   data = data_tmp;
   device = dstDevice;

   faust_cudaSetDevice(currentGPU);
}

template <typename faust_real>
void faust_cu_mat<faust_real>::setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_)
{
    resize(0,0,device_);
    dim1 = nbRow;
    dim2 = nbCol;
    isZeros = true;
    isIdentity = false;
}

template <typename faust_real>
void faust_cu_mat<faust_real>::setEyes(const faust_unsigned_int nbRow, const int device_)
{
    if(nbRow==device_)
        cerr << __FILE__ << ":" << __LINE__ << " : Warning - prototype is setEyes(const faust_unsigned_int nbRow, const int device_) and not setEyes(const int nbRow, const int nbCol). Please check whether syntax is correct" << endl;
    setZeros(nbRow, nbRow, device_);
    isIdentity = true;
    isZeros = false;
}

template <typename faust_real>
void faust_cu_mat<faust_real>::setEyes()
{
   if(dim1!=dim2)
      handleError(class_name,"setEyes() : GPU matrix is not square");
   setEyes(dim1, device);
}
 /// OPERATION BASIQUE ///

template <typename faust_real>
void faust_cu_mat<faust_real>::transpose(cublasHandle_t cublasHandle)
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

      faust_real* data_tmp;
      faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(faust_real));
      faust_cudaMemcpy(data_tmp, data, dim1*dim2*sizeof(faust_real), cudaMemcpyDeviceToDevice);

      faust_real alpha=1.0f;
      faust_real beta=0.0f;
      // According to the cublas doc, in-place mode if C = B, ldc = ldb and transb = CUBLAS_OP_N
      faust_cu_geam(cublasHandle,
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

template <typename faust_real>
void faust_cu_mat<faust_real>::init_from_transpose(const faust_cu_mat<faust_real>& cu_A, cublasHandle_t cublasHandle)
{
#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif

if(cu_A == (*this))
            handleError(class_name, "init_from_transpose(const faust_cu_mat<faust_real>&) : input GPU data is the same that the current GPU data. Try using faust_cu_mat<faust_real>::transpose instead, in this particular case.");

   if(cu_A.isZeros)
   {
      setZeros(cu_A.dim2, cu_A.dim1);
      #ifdef __COMPILE_TIMERS__
         t_transpose.stop();
      #endif
      return;
   }
   if(cu_A.isIdentity)
   {
      setEyes(cu_A.dim1);
      #ifdef __COMPILE_TIMERS__
         t_transpose.stop();
      #endif
      return;
   }
   resize(cu_A.dim2, cu_A.dim1);
   if (cu_A.dim1*cu_A.dim2 > 0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      faust_real* data_tmp = NULL;
      if(device == cu_A.device)
          data_tmp = cu_A.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, cu_A.dim1*cu_A.dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(data_tmp, device, cu_A.data, cu_A.device, cu_A.dim1*cu_A.dim2*sizeof(faust_real));
      }

      faust_real alpha=1.0f;
      faust_real beta=0.0f;
      faust_cu_geam(cublasHandle,
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
t_transpose.stop();
#endif

}

#ifdef __COMPILE_SPMAT__
template <typename faust_real>
void faust_cu_mat<faust_real>::init_from_transpose(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle)
{
#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif

    if(cu_S.getNbRow()*cu_S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.getNbRow(), cu_S.getNbCol());
       return;
    }

    if(cu_S.getNonZeros() == 0)
    {
       // display warning as matrix is zeros
       setZeros(cu_S.getNbRow(), cu_S.getNbCol());
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
            faust_cu_csc2dense(cusparseHandle,
               cu_S.getNbCol(), cu_S.getNbRow(), cu_S.getDescr(),  
               cu_S.getValues(), cu_S.getColInd(), cu_S.getRowPtr(),
               data, cu_S.getNbCol());
         else
         {
            int *rowptr_tmp, *colind_tmp;
            faust_real *values_tmp;
            faust_cudaMalloc((void**)&rowptr_tmp, (cu_S.getNbRow()+1)*sizeof(int));
            faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
            faust_cudaMalloc((void**)&values_tmp, cu_S.getNonZeros()*sizeof(faust_real));
            faust_cudaMemcpyPeer(rowptr_tmp, device, cu_S.getRowPtr(), cu_S.getDevice(), (cu_S.getNbRow()+1)*sizeof(int));
            faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
            faust_cudaMemcpyPeer(values_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(faust_real));
            faust_cu_csc2dense(cusparseHandle,
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
t_transpose.stop();
#endif

}
#endif



/*
 // former definition before multiply become an inline member
template <typename faust_real>
void faust_cu_mat<faust_real>::multiplyRight(const faust_cu_mat<faust_real>& cu_B, cublasHandle_t cublasHandle)
{

#ifdef __COMPILE_TIMERS__
t_mult_right.start();
#endif


   if (dim2 != cu_B.dim1)
   {
      handleError(class_name, "multiplyRight : dimension conflict between matrix");
   }

   if(cu_B.isIdentity)
   {
      #ifdef __COMPILE_TIMERS__
         t_mult_right.stop();
      #endif
      return;
   }

   if(isZeros || cu_B.isZeros)
   {
      //std::cout<<"zero"<<std::endl;
      setZeros(dim1, cu_B.dim2);
      isZeros = true;
      isIdentity = false;
      #ifdef __COMPILE_TIMERS__
         t_mult_right.stop();
      #endif
      return;
   }

   if(isIdentity)
   {
      this->operator=(cu_B);
      #ifdef __COMPILE_TIMERS__
         t_mult_right.stop();
      #endif
      return;
   }

     int currentGPU;
     faust_cudaGetDevice(&currentGPU);
     faust_cudaSetDevice(device);

      faust_real alpha = (faust_real)1.0;
      faust_real beta = (faust_real)0.0;

      faust_cu_mat<faust_real> cu_A(*this, device);
      faust_real* cu_B_copy;
      if (cu_B.device != device)
      {
         faust_cudaMalloc((void**)&cu_B_copy, cu_B.dim1*cu_B.dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(cu_B_copy, device, cu_B.getData(), cu_B.device, cu_B.dim1*cu_B.dim2*sizeof(faust_unsigned_int));
      }
      else
         cu_B_copy = cu_B.getData();

      resize(cu_A.dim1, cu_B.dim2);
      	faust_cu_gemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, dim1, cu_B.dim2, cu_A.dim2, &alpha, cu_A.data, cu_A.dim1, cu_B_copy, cu_B.dim1, &beta, getData(), dim1 );
         //cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim1, dim2, cu_A.dim2, 1.0f, cu_A.data, cu_A.dim1, A.getData(), A.dim1, 0.0f, getData(), dim1);
     if (cu_B.device != device)
     {
        faust_cudaFree(cu_B_copy);
        cu_B_copy = NULL;
     }
#ifdef __COMPILE_TIMERS__
t_mult_right.stop();
#endif

}
*/

/*
 // former definition before multiply become an inline member
template <typename faust_real>
void faust_cu_mat<faust_real>::multiplyLeft(faust_cu_mat<faust_real>& cu_A, cublasHandle_t cublasHandle)
{

#ifdef __COMPILE_TIMERS__
 t_mult_left.start();
#endif


   if (dim1 != cu_A.dim2)
      handleError(class_name, "multiplyLeft : dimension conflict between matrix");

   if(cu_A.isIdentity)
   {
      #ifdef __COMPILE_TIMERS__
         t_mult_left.stop();
      #endif
      return;
   }

   if(isZeros || cu_A.isZeros)
   {
      resize(cu_A.dim1,dim2);
      faust_real *const ptr_data_dst = getData();
      memset(ptr_data_dst, 0, sizeof(faust_real) * dim1*dim2);
      isZeros = true;
      #ifdef __COMPILE_TIMERS__
         t_mult_left.stop();
      #endif
      return;
   }

   if(isIdentity)
   {
      this->operator=(cu_A);
      #ifdef __COMPILE_TIMERS__
         t_mult_left.stop();
      #endif
      return;
   }


     int currentGPU;
     faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real alpha = (faust_real)1.0;
      faust_real beta = (faust_real)0.0;

      faust_cu_mat<faust_real> cu_B(*this, device);
      faust_real* cu_A_data;
      if (cu_A.device != device)
      {
         faust_cudaMalloc((void**)&cu_A_data, cu_A.dim1*cu_A.dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(cu_A_data, device, cu_A.getData(), cu_A.device, cu_A.dim1*cu_A.dim2*sizeof(faust_real));
      }
      else
         cu_A_data = cu_A.getData();

      resize(cu_A.dim1, cu_B.dim2);
      faust_cu_gemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, dim1, dim2, cu_A.dim2, &alpha, cu_A.data, cu_A.dim1, cu_B.data, cu_B.dim1, &beta, getData(), dim1 );
     if (cu_A.device != device)
     {
        faust_cudaFree(cu_A_data);
        cu_A_data = NULL;
     }

#ifdef __COMPILE_TIMERS__
 t_mult_left.stop();
#endif

}
*/

template <typename faust_real>
faust_real faust_cu_mat<faust_real>::max() const
{
   if(isZeros)
      return (faust_real)0.0;
   else if (isIdentity)
      return (faust_real)1.0;
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val =  faust_cu_max(data, dim1*dim2);
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
faust_real faust_cu_mat<faust_real>::min() const
{
   if(isZeros || isIdentity)
      return (faust_real)0.0;
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val =  faust_cu_min(data, dim1*dim2);
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
void faust_cu_mat<faust_real>::abs()
{
   if(isZeros || isIdentity)
      return ;
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      kernel_abs(data, dim1*dim2);
      faust_cudaSetDevice(currentGPU);
      return;
   }
   else
   {
      // display warning as matrix is empty
      return;
   }
}

template <typename faust_real>
faust_real faust_cu_mat<faust_real>::norm() const
{
   if(isZeros)
      return (faust_real)0.0;
   else if(isIdentity)
      return (faust_real)sqrt(dim1);
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);
      faust_real val =  faust_cu_norm(data, dim1*dim2);
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
faust_real faust_cu_mat<faust_real>::trace() const
{
   if(dim1 != dim2)
      handleError(class_name, "norm : matrix must be square");

   if (isZeros)
      return (faust_real)0.0;
   else if (isIdentity)
      return (faust_real)dim1;
   else if (dim1*dim2>0)
   {
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      faust_cudaSetDevice(device);

      faust_cu_vec<faust_real> cu_diag(dim1);
      kernel_get_diag(cu_diag.getData(), data, dim1);
      faust_cudaSetDevice(currentGPU);

      return faust_cu_sum(cu_diag.getData(), dim1);
   }
   else
   {
      // display warning as matrix is empty
      return (faust_real)0.0;
   }
}


/*template <typename faust_real>
 faust_real faust_cu_mat<faust_real>::spectralNorm() const
 {
   #ifdef __COMPILE_TIMERS__
         t_spectral_norm.start();
   #endif

    faust_real res=mat.operatorNorm();

   #ifdef __COMPILE_TIMERS__
      t_spectral_norm.stop();
   #endif
   return res;
 }*/


template <typename faust_real>
 faust_real faust_cu_mat<faust_real>::spectralNorm(const faust_unsigned_int nbr_iter_max,faust_real threshold, faust_int & flag, cublasHandle_t cublasHandle) const
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


   faust_cu_mat<faust_real> AtA;
   if (nb_row <= nb_col)
      gemm<faust_real>((*this),(*this),AtA,1.,0,'N','T',cublasHandle);
   else
      gemm<faust_real>((*this),(*this),AtA,1.,0,'T','N', cublasHandle);



   faust_real  res=std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag,cublasHandle));



   #ifdef __COMPILE_TIMERS__
      t_spectral_norm2.stop();
   #endif
   return res;

}


template <typename faust_real>
void faust_cu_mat<faust_real>::scalarMultiply(const faust_real lambda)
{
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.start();
#endif
   if(dim1*dim2 == 0)
   {
      // display warning as matrix is empty
      return;
   }

   if(isZeros){

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
         kernel_memset(data, (faust_real)0.0, dim1*dim2);
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

template <typename faust_real>
bool faust_cu_mat<faust_real>::operator==(const faust_cu_mat<faust_real>& cu_M)const
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
template <typename faust_real>
void faust_cu_mat<faust_real>::operator=(const faust_cu_mat<faust_real>& cu_M)
{
    if(cu_M.dim1*cu_M.dim2 == 0)
    {
       // display warning as matrix is empty
       resize(cu_M.dim1, cu_M.dim2);
       data = NULL ;
       device = cu_M.device;
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
}

template <typename faust_real>
void faust_cu_mat<faust_real>::operator=(const faust_mat<faust_real>& M)
{

    if(M.getNbRow()*M.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(M.getNbRow(), M.getNbCol());
       data = NULL ;
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
}

#ifdef __COMPILE_SPMAT__
template <typename faust_real>
void faust_cu_mat<faust_real>::operator=(const faust_spmat<faust_real>& S)
{
    if(S.getNbRow()*S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(S.getNbRow(), S.getNbCol());
       data = NULL ;
       return;
    }
    if(S.getNonZeros()==0)
    {
        setZeros(S.getNbRow(), S.getNbCol());
        return;
    }

    faust_mat<faust_real> M;
    M=S;

    int currentGPU;
    faust_cudaGetDevice(&currentGPU);
    faust_cudaSetDevice(device);

    faust_cudaMemcpy(data, M.getData(), dim1*dim2*sizeof(faust_real), cudaMemcpyHostToDevice);

    faust_cudaSetDevice(currentGPU);

    isZeros = false;
    isIdentity = false;
}

template <typename faust_real>
void faust_cu_mat<faust_real>::operator=(const faust_cu_spmat<faust_real>& S)
{handleError(class_name, "faust_cu_mat<faust_real>::operator=(const faust_cu_spmat&) is not defined. Use faust_cu_mat<faust_real>::init_from_cu_spmat(const faust_cu_spmat&, cusparseHandle_t) instead");}

template <typename faust_real>
void faust_cu_mat<faust_real>::init_from_cu_spmat(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle, const faust_real coeff /*=1.0*/)
{
    if(cu_S.getNbRow()*cu_S.getNbCol() == 0)
    {
       // display warning as matrix is empty
       resize(cu_S.getNbRow(), cu_S.getNbCol());
       return;
    }

    if(cu_S.getNonZeros() == 0)
    {
       // display warning as matrix is zeros
       setZeros(cu_S.getNbRow(), cu_S.getNbCol());
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
               faust_cu_csr2dense(cusparseHandle,
                  cu_S.getNbRow(), cu_S.getNbCol(), cu_S.getDescr(),  
                  cu_S.getValues(), cu_S.getRowPtr(), cu_S.getColInd(),
                  data, cu_S.getNbRow());
         }
         else
         {
            faust_real *values_tmp;
            faust_cudaMalloc((void**)&values_tmp, cu_S.getNonZeros()*sizeof(faust_real));
            faust_cudaMemcpyPeer(values_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(faust_real));
            if(coeff != 1.0)
               kernel_mult_const(values_tmp, coeff, cu_S.getNonZeros());
            
            if(device!=cu_S.getDevice())
            {
               int *rowptr_tmp, *colind_tmp;

               faust_cudaMalloc((void**)&rowptr_tmp, cu_S.getNonZeros()*sizeof(int));
               faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));

               faust_cudaMemcpyPeer(rowptr_tmp, device, cu_S.getRowPtr(), cu_S.getDevice(), (cu_S.getNbRow()+1)*sizeof(int));
               faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));


               faust_cu_csr2dense(cusparseHandle,
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
               faust_cu_csr2dense(cusparseHandle,
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
}

template <typename faust_real>
void faust_cu_mat<faust_real>::operator+=(const faust_cu_spmat<faust_real>& cu_S)
{handleError(class_name, "faust_cu_mat<faust_real>::operator+=(const faust_cu_spmat&) is not defined. Use faust_cu_mat<faust_real>::add(const faust_cu_spmat&, cusparseHandle_t) instead");}

template <typename faust_real>
void faust_cu_mat<faust_real>::add(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle)
{
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
      init_from_cu_spmat(cu_S, cusparseHandle);
   else if (cu_S.getNonZeros()==0){} // We don't do anything but it is written not to go to else
   else if(isIdentity)
   {
      init_from_cu_spmat(cu_S, cusparseHandle);
      kernel_add_diag_const(data, (faust_real)1.0, dim1);
   }
   else
   {

      faust_cudaSetDevice(cu_S.getDevice());

      int* rowind_tmp_src;
      faust_cudaMalloc((void**)&rowind_tmp_src, cu_S.getNonZeros()*sizeof(int));
      faust_cu_csr2coo(cusparseHandle, cu_S.getRowPtr(),
         cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaSetDevice(device);
       

      if(device == cu_S.getDevice())
      {
         faust_cu_csr2coo(cusparseHandle, cu_S.getRowPtr(),
            cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

         kernel_add_sparse2full(data, rowind_tmp_src, cu_S.getColInd(), cu_S.getValues(), cu_S.getNonZeros(), dim1);

      }
      else
      {
         int *rowind_tmp_dst, *colind_tmp;
         faust_real* data_tmp;

         faust_cudaMalloc((void**)&rowind_tmp_dst, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&data_tmp, cu_S.getNonZeros()*sizeof(faust_real));
         faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(rowind_tmp_dst, device, rowind_tmp_src, cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(data_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(faust_real));
         kernel_add_sparse2full(data, rowind_tmp_dst, colind_tmp, data_tmp, cu_S.getNonZeros(), dim1);
         faust_cudaFree(rowind_tmp_dst);
         faust_cudaFree(colind_tmp);
         faust_cudaFree(data_tmp);
      }
      faust_cudaFree(rowind_tmp_src);
   }

   faust_cudaSetDevice(currentGPU);  
}


template <typename faust_real>
void faust_cu_mat<faust_real>::operator-=(const faust_cu_spmat<faust_real>& cu_S)
{handleError(class_name, "faust_cu_mat<faust_real>::operator-=(const faust_cu_spmat&) is not defined. Use faust_cu_mat<faust_real>::sub(const faust_cu_spmat&, cusparseHandle_t) instead");}

template <typename faust_real>
void faust_cu_mat<faust_real>::sub(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle)
{
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
      init_from_cu_spmat(cu_S, cusparseHandle, -1.0);
   else if (cu_S.getNonZeros()==0){} // We don't do anything but it is written not to go to else
   else if(isIdentity)
   {
      init_from_cu_spmat(cu_S, cusparseHandle, -1.0);
      kernel_add_diag_const(data, (faust_real)1.0, dim1);
   }
   else
   {
       // if *this and cu_S are the same object, we can use kernel_add, so there is no need to know if these matrix are the same objet or not

      faust_cudaSetDevice(cu_S.getDevice());

      int* rowind_tmp_src;
      cudaMalloc((void**)&rowind_tmp_src, cu_S.getNonZeros()*sizeof(int));
      faust_cu_csr2coo(cusparseHandle, cu_S.getRowPtr(),
         cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

      faust_cudaSetDevice(device);
       

      if(device == cu_S.getDevice())
      {
         faust_cu_csr2coo(cusparseHandle, cu_S.getRowPtr(),
            cu_S.getNonZeros(), dim1, rowind_tmp_src, CUSPARSE_INDEX_BASE_ZERO);

         kernel_sub_sparse2full(data, rowind_tmp_src, cu_S.getColInd(), cu_S.getValues(), cu_S.getNonZeros(), dim1);

      }
      else
      {
         int *rowind_tmp_dst, *colind_tmp;
         faust_real* data_tmp;

         faust_cudaMalloc((void**)&rowind_tmp_dst, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&colind_tmp, cu_S.getNonZeros()*sizeof(int));
         faust_cudaMalloc((void**)&data_tmp, cu_S.getNonZeros()*sizeof(faust_real));
         faust_cudaMemcpyPeer(colind_tmp, device, cu_S.getColInd(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(rowind_tmp_dst, device, rowind_tmp_src, cu_S.getDevice(), cu_S.getNonZeros()*sizeof(int));
         faust_cudaMemcpyPeer(data_tmp, device, cu_S.getValues(), cu_S.getDevice(), cu_S.getNonZeros()*sizeof(faust_real));
         kernel_sub_sparse2full(data, rowind_tmp_dst, colind_tmp, data_tmp, cu_S.getNonZeros(), dim1);
         faust_cudaFree(rowind_tmp_dst);
         faust_cudaFree(colind_tmp);
         faust_cudaFree(data_tmp);
      }
      faust_cudaFree(rowind_tmp_src);
   }

   faust_cudaSetDevice(currentGPU);  
}
#endif //#ifdef __COMPILE_SPMAT__


template <typename faust_real>
void faust_cu_mat<faust_real>::add(const faust_cu_mat<faust_real>& cu_M)
{
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
      kernel_memset(data, (faust_real)0.0, dim1*dim2);
      kernel_add_diag_const(data, (faust_real)2.0, dim1);
   }
   else if(isIdentity)
   {
      this->operator=(cu_M);
      kernel_add_diag_const(data, (faust_real)1.0, dim1);
   }
   else if(cu_M.isIdentity)
      kernel_add_diag_const(data, (faust_real)1.0, dim1);

   else
   {
       // if *this and cu_M are the same object, we can use kernel_add, so there is no need to know if these matrix are the same objet or not
      if(device == cu_M.device)
         kernel_add(data, cu_M.data, dim1*dim2);
      else
      {
         faust_real* data_tmp;
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(faust_real));
         kernel_add(data, data_tmp, dim1*dim2);
         faust_cudaFree(data_tmp);
      }
   }

   faust_cudaSetDevice(currentGPU);

}



template <typename faust_real>
void faust_cu_mat<faust_real>::operator-=(const faust_cu_mat<faust_real>& cu_M)
{
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
      kernel_add_diag_const(data, (faust_real)-1.0, dim1);
   else
   {
      faust_real* data_tmp = NULL;
      if(device == cu_M.device)
         data_tmp = cu_M.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(faust_real));
      }

      if (isZeros)
      {
         resize(dim1, dim2);
         kernel_memset(data, (faust_real)0.0, dim1*dim2);
         kernel_sub(data, data_tmp, dim1*dim2);
      }
      else if(isIdentity)
      {
         resize(dim1, dim2);
         kernel_memset(data, (faust_real)0.0, dim1*dim2);
         kernel_sub(data, data_tmp, dim1*dim2);
         kernel_add_diag_const(data, (faust_real)1.0, dim1);
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

}



template <typename faust_real>
void faust_cu_mat<faust_real>::scalarMultiply(const faust_cu_mat<faust_real>& cu_M)
{
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
      faust_cu_mat<faust_real> cu_A_tmp(*this);
      kernel_memset(data, (faust_real)0.0, dim1*dim2);
      kernel_copy_diag(data, cu_A_tmp.data, dim1);
   }
   else
   {
      faust_real* data_tmp = NULL;
      if(device == cu_M.device)
         data_tmp = cu_M.data;
      else
      {
         faust_cudaMalloc((void**)&data_tmp, dim1*dim2*sizeof(faust_real));
         faust_cudaMemcpyPeer(data_tmp, device, cu_M.data, cu_M.device, dim1*dim2*sizeof(faust_real));
      }

      if(isIdentity)
      {
         kernel_memset(data, (faust_real)0.0, dim1*dim2);
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
}

template <typename faust_real>
void faust_cu_mat<faust_real>::Display()const
{
   faust_mat<faust_real> M;
   faust_cu2faust(M,*this);
   M.Display();
}

template <typename faust_real>
void faust_cu_mat<faust_real>::print_file(const char* filename)const
{
   faust_mat<faust_real> M;
   faust_cu2faust(M,*this);
   M.print_file(filename);
}



#ifdef __COMPILE_TIMERS__
faust_timer faust_cu_mat<faust_real>::t_constr;
faust_timer faust_cu_mat<faust_real>::t_get_coeff;
faust_timer faust_cu_mat<faust_real>::t_get_coeffs;
faust_timer faust_cu_mat<faust_real>::t_set_coeff;
faust_timer faust_cu_mat<faust_real>::t_set_coeffs;
faust_timer faust_cu_mat<faust_real>::t_set_coeffs2;
faust_timer faust_cu_mat<faust_real>::t_resize;
faust_timer faust_cu_mat<faust_real>::t_check_dim;
faust_timer faust_cu_mat<faust_real>::t_max;
faust_timer faust_cu_mat<faust_real>::t_transpose;
faust_timer faust_cu_mat<faust_real>::t_mult_right;
faust_timer faust_cu_mat<faust_real>::t_mult_left;
faust_timer faust_cu_mat<faust_real>::t_scalar_multiply;
faust_timer faust_cu_mat<faust_real>::t_add;
faust_timer faust_cu_mat<faust_real>::t_sub;
faust_timer faust_cu_mat<faust_real>::t_print_file;
faust_timer faust_cu_mat<faust_real>::t_multiply;
faust_timer faust_cu_mat<faust_real>::t_gemm;
faust_timer faust_cu_mat<faust_real>::t_add_ext;

faust_timer faust_cu_mat<faust_real>::t_spectral_norm;
faust_timer faust_cu_mat<faust_real>::t_spectral_norm2;
faust_timer faust_cu_mat<faust_real>::t_power_iteration;

template <typename faust_real>
void faust_cu_mat<faust_real>::print_timers()const
{
   cout << "timers in faust_cu_mat :" << endl;
   cout << "t_constr          = " << t_constr.get_time()          << " s for "<< t_constr.get_nb_call()          << " calls" << endl;
   cout << "t_get_coeff       = " << t_get_coeff.get_time()       << " s for "<< t_get_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_get_coeffs      = " << t_get_coeffs.get_time()      << " s for "<< t_get_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeff       = " << t_set_coeff.get_time()       << " s for "<< t_set_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_set_coeffs      = " << t_set_coeffs.get_time()      << " s for "<< t_set_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeffs2     = " << t_set_coeffs2.get_time()     << " s for "<< t_set_coeffs2.get_nb_call()     << " calls" << endl;
   cout << "t_resize          = " << t_resize.get_time()          << " s for "<< t_resize.get_nb_call()          << " calls" << endl;
   cout << "t_check_dim       = " << t_check_dim.get_time()       << " s for "<< t_check_dim.get_nb_call()       << " calls" << endl;
   cout << "t_max             = " << t_max.get_time()             << " s for "<< t_max.get_nb_call()             << " calls" << endl;
   cout << "t_transpose       = " << t_transpose.get_time()       << " s for "<< t_transpose.get_nb_call()       << " calls" << endl;
   cout << "t_mult_right      = " << t_mult_right.get_time()      << " s for "<< t_mult_right.get_nb_call()      << " calls" << endl;
   cout << "t_mult_left       = " << t_mult_left.get_time()       << " s for "<< t_mult_left.get_nb_call()       << " calls" << endl;
  cout << "t_scalar_multiply = " << t_scalar_multiply.get_time() << " s for "<< t_scalar_multiply.get_nb_call() << " calls" << endl;
   cout << "t_add             = " << t_add.get_time()             << " s for "<< t_add.get_nb_call()             << " calls" << endl;
   cout << "t_sub             = " << t_sub.get_time()             << " s for "<< t_sub.get_nb_call()             << " calls" << endl;
cout << "t_print_file      = " << t_print_file.get_time()      << " s for "<< t_print_file.get_nb_call()      << " calls" << endl<<endl;


   cout << "timers in faust_cu_mat / LinearAlgebra :" << endl;
   cout << "t_multiply        = " << t_multiply.get_time()        << " s for "<< t_multiply.get_nb_call()        << " calls" << endl;
   cout << "t_gemm            = " << t_gemm.get_time()            << " s for "<< t_gemm.get_nb_call()            << " calls" << endl;
   cout << "t_add_ext         = " << t_add_ext.get_time()         << " s for "<< t_add_ext.get_nb_call()         << " calls" << endl<<endl<<endl;
}
#endif

#endif
