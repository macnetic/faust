/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_CUDA_H__
#define __FAUST_CUDA_H__

#include <iostream>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "cublas_v2.h"
#include "cusparse.h"


inline void faust_cu_csr2coo(cusparseHandle_t handle, const int *csrRowPtr,
   int nnz, int m, int *cooRowInd, cusparseIndexBase_t idxBase);

template <typename faust_real>
void faust_cu_dot(cublasHandle_t handle, int n,
   const faust_real* x, int incx, const faust_real* y, int incy, faust_real* result);


template <typename faust_real>
void faust_cu_gemv(cublasHandle_t handle, cublasOperation_t trans,
   int m, int n, const faust_real* alpha, const faust_real* A, int lda,
   const faust_real* x, int incx, const faust_real* beta, faust_real* y, int incy);

template <typename faust_real>
void faust_cu_gemm(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n, int k,
   const faust_real* alpha, const faust_real* A, int lda,
   const faust_real* B, int ldb,
   const faust_real* beta, faust_real*C , int ldc);

template <typename faust_real>
void faust_cu_gemm(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n,
   const faust_real* alpha, const faust_real* A, int lda,
   const faust_real* beta, const faust_real* B, int ldb,
   faust_real* C, int ldc);

template <typename faust_real>
void faust_cu_geam(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n,
   const faust_real* alpha, const faust_real* A, int lda,
   const faust_real* beta, const faust_real* B, int ldb,
   faust_real* C, int ldc);



template <typename faust_real>
void faust_cu_csrmv(cusparseHandle_t handle, cusparseOperation_t transA,
   int m, int n, int nnz, const faust_real* alpha, const cusparseMatDescr_t descrA,
   const faust_real* csrValA,  const int* csrRowPtrA, const int *csrColIndA,
   const faust_real* x, const faust_real* beta, faust_real* y);

template <typename faust_real>
void faust_cu_csrmm(cusparseHandle_t handle, 
   cusparseOperation_t transA, cusparseOperation_t transB,    
   int m, int n, int k, int nnz,       
   const faust_real* alpha, const cusparseMatDescr_t descrA,
   const faust_real* csrValA, const int* csrRowPtrA, const int* csrColIndA,
   const faust_real* B, int ldb, const faust_real* beta, faust_real* C, int ldc);

template <typename faust_real> void faust_cu_nnz(cusparseHandle_t handle, cusparseDirection_t dirA, 
   int m, int n, const cusparseMatDescr_t descrA, const faust_real* A, 
   int lda, int* nnzPerRowColumn, int* nnzTotalDevHostPtr);

template <typename faust_real>
void faust_cu_dense2csr(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const faust_real* A, int lda, const int* nnzPerRow,
   faust_real* csrValA, int* csrRowPtrA, int* csrColIndA);

template <typename faust_real>
void faust_cu_csr2dense(cusparseHandle_t handle,
 int m, int n, const cusparseMatDescr_t descrA,  
 const faust_real* csrValA, const int* csrRowPtrA, 
 const int* csrColIndA, faust_real* A, int lda);

template <typename faust_real>
void faust_cu_csc2dense(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const faust_real* cscValA, const int* cscRowIndA, 
   const int* cscColPtrA, faust_real* A, int lda);

template <typename faust_real>
void faust_cu_csr2csc(cusparseHandle_t handle, int m, int n, int nnz,
   const faust_real* csrVal, const int* csrRowPtr, const int* csrColInd, 
   faust_real* cscVal, int* cscRowInd, int* cscColPtr,
   cusparseAction_t copyValues, cusparseIndexBase_t idxBase);


/*
#define cutilSC(err,commandName)           __cudaSC      (err, __FILE__, __LINE__,commandName)
inline void __cudaSC(cudaError err, const char *file, const int line, const char* cmdName) {
  if (rudaSuccess != err) {
    std::cerr<<"Error in "<<file<<":"<<line<<" : cmdName failed : "<<cudaGetErrorString(err)<<std::endl;
    exit(EXIT_FAILURE);
  }
}*/

inline const char* cublasGetErrorString(cublasStatus_t cublasStat)
{
    switch(cublasStat)
    {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE"; 
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH"; 
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED"; 
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR"; 
        case CUBLAS_STATUS_NOT_SUPPORTED: return "CUBLAS_STATUS_NOT_SUPPORTED"; 
        case CUBLAS_STATUS_LICENSE_ERROR: return "CUBLAS_STATUS_LICENSE_ERROR"; 
    }
    return "unknown error";
} 


inline const char* cusparseGetErrorString(cusparseStatus_t cusparseStat)
{
    switch(cusparseStat)
    {
        case CUSPARSE_STATUS_SUCCESS: return "CUSPARSE_STATUS_SUCCESS";
        case CUSPARSE_STATUS_NOT_INITIALIZED: return "CUSPARSE_STATUS_NOT_INITIALIZED";
        case CUSPARSE_STATUS_ALLOC_FAILED: return "CUSPARSE_STATUS_ALLOC_FAILED";
        case CUSPARSE_STATUS_INVALID_VALUE: return "CUSPARSE_STATUS_INVALID_VALUE";
        case CUSPARSE_STATUS_ARCH_MISMATCH: return "CUSPARSE_STATUS_ARCH_MISMATCH";
        case CUSPARSE_STATUS_MAPPING_ERROR: return "CUSPARSE_STATUS_MAPPING_ERROR";
        case CUSPARSE_STATUS_EXECUTION_FAILED: return "CUSPARSE_STATUS_EXECUTION_FAILED";
        case CUSPARSE_STATUS_INTERNAL_ERROR: return "CUSPARSE_STATUS_INTERNAL_ERROR";
        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED: return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    }    
    return "unknown error";
}


#define faust_cudaSafe(call, commandName) do {                                \
    cudaError err = call;                                                    \
    if( err != cudaSuccess) {                                                \
        std::cerr << __FILE__<<":"<<__LINE__<< " : Error : " << commandName " failed : " << cudaGetErrorString( err) << std::endl ;              \
        exit(err);                                                  \
    } } while (0)

#define faust_cublasSafe(call, commandName) do {                                \
    cublasStatus_t cublasStat = call;                                                    \
    if( cublasStat != CUBLAS_STATUS_SUCCESS) {                                                \
        std::cerr << __FILE__<<":"<<__LINE__<< " : Error : " << commandName " failed : " << cublasGetErrorString(cublasStat) << std::endl ;              \
        exit(cublasStat);                                                  \
    } } while (0)

#define faust_cusparseSafe(call, commandName) do {                                \
    cusparseStatus_t cusparseStat = call;                                                    \
    if( cusparseStat != CUSPARSE_STATUS_SUCCESS) {                                                \
        std::cerr << __FILE__<<":"<<__LINE__<< " : Error : " << commandName " failed : " << cusparseGetErrorString(cusparseStat) << std::endl ;              \
        exit(cusparseStat);                                                  \
    } } while (0)

#define faust_kernelSafe() do {                                \
    cudaError_t err = cudaGetLastError();                                                    \
    if( err != cudaSuccess) {                                                \
        std::cerr << __FILE__<<":"<<__LINE__<< " : Error : kernel failed : " << cudaGetErrorString(cudaSuccess) << std::endl ;              \
        exit(err);                                                  \
    } } while (0)




#define faust_cudaGetDevice(arg) do {                                \
 faust_cudaSafe(cudaGetDevice(arg),"cudaGetDevice");  \
     } while (0)

#define faust_cudaSetDevice(arg) do {                                \
 faust_cudaSafe(cudaSetDevice(arg),"cudaSetDevice");  \
     } while (0)

#define faust_cudaFree(arg) do {                                \
 faust_cudaSafe(cudaFree(arg),"cudaFree");  \
     } while (0)

#define faust_cudaMalloc(ptr_address, byte_size) do {                                \
 faust_cudaSafe(cudaMalloc(ptr_address, byte_size),"cudaMalloc");  \
     } while (0)

#define faust_cudaMemcpy(data_dst, data_src, byte_size, direction) do {                                \
 faust_cudaSafe(cudaMemcpy(data_dst, data_src, byte_size, direction),"cudaMemcpy");  \
     } while (0)

#define faust_cudaMemcpyPeer(data_dst, dev_dst, data_src, dev_src, byte_size) do {                                \
 faust_cudaSafe(cudaMemcpyPeer(data_dst, dev_dst, data_src, dev_src, byte_size),"cudaMemcpyPeer");  \
     } while (0)

#define faust_cudaMemcpyAsync(data_dst, data_src, byte_size, direction, stream) do {                                \
 faust_cudaSafe(cudaMemcpyAsync(data_dst, data_src, byte_size, direction, stream),"cudaMemcpyAsync");  \
     } while (0)

#define faust_cudaMemcpyPeerAsync(data_dst, dev_dst, data_src, dev_src, byte_size, stream) do {  \
 faust_cudaSafe(cudaMemcpyPeerAsync(data_dst, dev_dst, data_src, dev_src, byte_size, stream),"cudaMemcpyPeerAsync");  \
     } while (0)

#define faust_cudaEventRecord(event, stream) do {  \
 faust_cudaSafe(cudaEventRecord(event, stream),"cudaEventRecord");  \
     } while (0)

#define faust_cudaEventSynchronize(event) do {  \
 faust_cudaSafe(cudaEventSynchronize(event),"cudaEventSynchronize");  \
     } while (0)

#define faust_cudaEventElapsedTime(ms, start, end) do {  \
 faust_cudaSafe(cudaEventElapsedTime(ms, start, end),"cudaEventElapsedTime");  \
     } while (0)

#define faust_cudaEventCreate(event) do {  \
 faust_cudaSafe(cudaEventCreate(event),"cudaEventCreate");  \
     } while (0)





#define faust_cublasCreate(handle) do {                  \
 faust_cublasSafe(cublasCreate(handle),"cublasCreate");  \
     } while (0)

#define faust_cublasDestroy(handle) do {                   \
 faust_cublasSafe(cublasDestroy(handle),"cublasDestroy");  \
     } while (0)

#define faust_cublasSdot(handle,n,x,incx,y,incy,result) do {                   \
 faust_cublasSafe(cublasSdot(handle,n,x,incx,y,incy,result),"cublasSdot");  \
     } while (0)
#define faust_cublasDdot(handle,n,x,incx,y,incy,result) do {                   \
 faust_cublasSafe(cublasDdot(handle,n,x,incx,y,incy,result),"cublasDdot");  \
     } while (0)

#define faust_cublasSgemv(handle,opA,m,n,alpha,A_data,lda,x_data,incx,beta,y_data,incy) do {                   \
 faust_cublasSafe(cublasSgemv(handle,opA,m,n,alpha,A_data,lda,x_data,incx,beta,y_data,incy),"cublasSgemv");  \
     } while (0)
#define faust_cublasDgemv(handle,opA,m,n,alpha,A_data,lda,x_data,incx,beta,y_data,incy) do {                   \
 faust_cublasSafe(cublasDgemv(handle,opA,m,n,alpha,A_data,lda,x_data,incx,beta,y_data,incy),"cublasDgemv");  \
     } while (0)

#define faust_cublasSgemm(handle,opA,opB,m,n,k,alpha,A_data,lda,B_data,ldB,beta,C_data,ldc) do {                   \
 faust_cublasSafe(cublasSgemm(handle,opA,opB,m,n,k,alpha,A_data,lda,B_data,ldB,beta,C_data,ldc),"cublasSgemm");  \
     } while (0)
#define faust_cublasDgemm(handle,opA,opB,m,n,k,alpha,A_data,lda,B_data,ldB,beta,C_data,ldc) do {                   \
 faust_cublasSafe(cublasDgemm(handle,opA,opB,m,n,k,alpha,A_data,lda,B_data,ldB,beta,C_data,ldc),"cublasDgemm");  \
     } while (0)

#define faust_cublasSgeam(handle,opA,opB,m,n,alpha,A_data,lda,beta,B_data,ldb,C_data,ldc) do {                   \
 faust_cublasSafe(cublasSgeam(handle,opA,opB,m,n,alpha,A_data,lda,beta,B_data,ldb,C_data,ldc),"cublasSgeam");  \
     } while (0)
#define faust_cublasDgeam(handle,opA,opB,m,n,alpha,A_data,lda,beta,B_data,ldb,C_data,ldc) do {                   \
 faust_cublasSafe(cublasDgeam(handle,opA,opB,m,n,alpha,A_data,lda,beta,B_data,ldb,C_data,ldc),"cublasDgeam");  \
     } while (0)




#define faust_cusparseCreate(handle) do {                  \
 faust_cusparseSafe(cusparseCreate(handle),"cusparseCreate");  \
     } while (0)

#define faust_cusparseCreateMatDescr(descrA) do {                  \
 faust_cusparseSafe(cusparseCreateMatDescr(descrA),"cusparseCreateMatDescr");  \
     } while (0)

#define faust_cusparseDestroyMatDescr(descrA) do {                  \
 faust_cusparseSafe(cusparseDestroyMatDescr(descrA),"cusparseDestroyMatDescr");  \
     } while (0)

#define faust_cusparseDestroy(handle) do {                  \
 faust_cusparseSafe(cusparseDestroy(handle),"cusparseDestroy");  \
     } while (0)

#define faust_cusparseScsrmv(handle,opA,m,n,nnz,alpha,A_descr,A_Values,A_RowPtr,A_ColInd,x_data,beta,y_data) do {                  \
 faust_cusparseSafe(cusparseScsrmv(handle,opA,m,n,nnz,alpha,A_descr,A_Values,A_RowPtr,A_ColInd,x_data,beta,y_data),"cusparseScsrmv");  \
     } while (0)
#define faust_cusparseDcsrmv(handle,opA,m,n,nnz,alpha,A_descr,A_Values,A_RowPtr,A_ColInd,x_data,beta,y_data) do {                  \
 faust_cusparseSafe(cusparseDcsrmv(handle,opA,m,n,nnz,alpha,A_descr,A_Values,A_RowPtr,A_ColInd,x_data,beta,y_data),"cusparseDcsrmv");  \
     } while (0)

#define faust_cusparseScsrmm2(handle,opA,opB,m,n,k,nnz,alpha,A_Descr,A_Values,A_RowPtr,A_ColInd,B_data,ldb,beta,C_data,ldc) do {                  \
 faust_cusparseSafe(cusparseScsrmm2(handle,opA,opB,m,n,k,nnz,alpha,A_Descr,A_Values,A_RowPtr,A_ColInd,B_data,ldb,beta,C_data,ldc),"cusparseScsrmm2");  \
     } while (0)
#define faust_cusparseDcsrmm2(handle,opA,opB,m,n,k,nnz,alpha,A_Descr,A_Values,A_RowPtr,A_ColInd,B_data,ldb,beta,C_data,ldc) do {                  \
 faust_cusparseSafe(cusparseDcsrmm2(handle,opA,opB,m,n,k,nnz,alpha,A_Descr,A_Values,A_RowPtr,A_ColInd,B_data,ldb,beta,C_data,ldc),"cusparseDcsrmm2");  \
     } while (0)

#define faust_cusparseSnnz(handle, dirA, m, n, descrA, A_data, lda, nnzPerRowCol, nnzTotal) do {                  \
 faust_cusparseSafe(cusparseSnnz(handle, dirA, m, n, descrA, A_data, lda, nnzPerRowCol, nnzTotal),"cusparseSnnz");  \
     } while (0)
#define faust_cusparseDnnz(handle, dirA, m, n, descrA, A_data, lda, nnzPerRowCol, nnzTotal) do {                  \
 faust_cusparseSafe(cusparseDnnz(handle, dirA, m, n, descrA, A_data, lda, nnzPerRowCol, nnzTotal),"cusparseDnnz");  \
     } while (0)

#define faust_cusparseSdense2csr(handle, m, n, descrA, A_data, lda, nnzPerRow, csrValues, csrRowPtr, csrColInd) do {                  \
 faust_cusparseSafe(cusparseSdense2csr(handle, m, n, descrA, A_data, lda, nnzPerRow, csrValues, csrRowPtr, csrColInd),"cusparseSdense2csr");  \
     } while (0)
#define faust_cusparseDdense2csr(handle, m, n, descrA, A_data, lda, nnzPerRow, csrValues, csrRowPtr, csrColInd) do {                  \
 faust_cusparseSafe(cusparseDdense2csr(handle, m, n, descrA, A_data, lda, nnzPerRow, csrValues, csrRowPtr, csrColInd),"cusparseDdense2csr");  \
     } while (0)

#define faust_cusparseScsr2dense(handle, m, n, descrA, csrValues, csrRowPtr, csrColInd, A_data, lda) do {                  \
 faust_cusparseSafe(cusparseScsr2dense(handle, m, n, descrA, csrValues, csrRowPtr, csrColInd, A_data, lda),"cusparseScsr2dense");  \
     } while (0)
#define faust_cusparseDcsr2dense(handle, m, n, descrA, csrValues, csrRowPtr, csrColInd, A_data, lda) do {                  \
 faust_cusparseSafe(cusparseDcsr2dense(handle, m, n, descrA, csrValues, csrRowPtr, csrColInd, A_data, lda),"cusparseScsr2dense");  \
     } while (0)

#define faust_cusparseScsc2dense(handle, m, n, descrA, cscValues, cscRowInd, cscColPtr, A_data, lda) do {                  \
 faust_cusparseSafe(cusparseScsc2dense(handle, m, n, descrA, cscValues, cscRowInd, cscColPtr, A_data, lda),"cusparseScsc2dense");  \
     } while (0)
#define faust_cusparseDcsc2dense(handle, m, n, descrA, cscValues, cscRowInd, cscColPtr, A_data, lda) do {                  \
 faust_cusparseSafe(cusparseDcsc2dense(handle, m, n, descrA, cscValues, cscRowInd, cscColPtr, A_data, lda),"cusparseDcsc2dense");  \
     } while (0)

#define faust_cusparseXcsr2coo(handle, csrRowPtr, nnz, m, cooRowInd, idxBase) do {                  \
 faust_cusparseSafe(cusparseXcsr2coo(handle, csrRowPtr, nnz, m, cooRowInd, idxBase),"cusparseXcsr2coo");  \
     } while (0)


#define faust_cusparseScsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase)  do {    \
 faust_cusparseSafe(cusparseScsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase),"cusparseScsr2csc");  \
     } while (0)
#define faust_cusparseDcsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase) do {    \
 faust_cusparseSafe(cusparseDcsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase),"cusparseDcsr2csc");  \
     } while (0)


#include "faust_cuda.hpp"

#endif

