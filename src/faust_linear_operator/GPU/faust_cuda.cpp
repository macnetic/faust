#include "faust_cuda.h"

template<>
void faust_cu_dot<float>(cublasHandle_t handle, int n,
   const float* x, int incx, const float* y, int incy, float* result)
{faust_cublasSdot(handle,n,x,incx,y,incy,result);}

template<>
void faust_cu_dot<double>(cublasHandle_t handle, int n,
   const double* x, int incx, const double* y, int incy, double* result)
{faust_cublasDdot(handle,n,x,incx,y,incy,result);}


template<>
void faust_cu_gemv<float>(cublasHandle_t handle, cublasOperation_t trans,
   int m, int n, const float* alpha, const float* A, int lda,
   const float* x, int incx, const float* beta, float* y, int incy)
{faust_cublasSgemv(handle,trans,m,n,alpha,A,lda,x,incx,beta,y,incy);}

template<>
void faust_cu_gemv<double>(cublasHandle_t handle, cublasOperation_t trans,
   int m, int n, const double* alpha, const double* A, int lda,
   const double* x, int incx, const double* beta, double* y, int incy)
{faust_cublasDgemv(handle,trans,m,n,alpha,A,lda,x,incx,beta,y,incy);}


template<>
void faust_cu_gemm<float>(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n, int k,
   const float* alpha, const float* A, int lda,
   const float* B, int ldb,
   const float* beta, float*C , int ldc)
{faust_cublasSgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);}

template<>
void faust_cu_gemm<double>(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n, int k,
   const double* alpha, const double* A, int lda,
   const double* B, int ldb,
   const double* beta, double*C , int ldc)
{faust_cublasDgemm(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);}


template<>
void faust_cu_geam<float>(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n,
   const float* alpha, const float* A, int lda,
   const float* beta, const float* B, int ldb,
   float* C, int ldc)
{faust_cublasSgeam(handle,transa,transb,m,n,alpha,A,lda,beta,B,ldb,C,ldc);}

template<>
void faust_cu_geam<double>(cublasHandle_t handle,
   cublasOperation_t transa, cublasOperation_t transb,
   int m, int n,
   const double* alpha, const double* A, int lda,
   const double* beta, const double* B, int ldb,
   double* C, int ldc)
{faust_cublasDgeam(handle,transa,transb,m,n,alpha,A,lda,beta,B,ldb,C,ldc);}




template<>
void faust_cu_csrmv<float>(cusparseHandle_t handle, cusparseOperation_t transA,
   int m, int n, int nnz, const float* alpha, const cusparseMatDescr_t descrA,
   const float* csrValA,  const int* csrRowPtrA, const int *csrColIndA,
   const float* x, const float* beta, float* y)
{faust_cusparseScsrmv(handle,transA,m,n,nnz,alpha,descrA,csrValA,csrRowPtrA,csrColIndA,x,beta,y);}

template<>
void faust_cu_csrmv<double>(cusparseHandle_t handle, cusparseOperation_t transA,
   int m, int n, int nnz, const double* alpha, const cusparseMatDescr_t descrA,
   const double* csrValA,  const int* csrRowPtrA, const int *csrColIndA,
   const double* x, const double* beta, double* y)
{faust_cusparseDcsrmv(handle,transA,m,n,nnz,alpha,descrA,csrValA,csrRowPtrA,csrColIndA,x,beta,y);}


template<>
void faust_cu_csrmm<float>(cusparseHandle_t handle, 
   cusparseOperation_t transA, cusparseOperation_t transB,    
   int m, int n, int k, int nnz,       
   const float* alpha, const cusparseMatDescr_t descrA,
   const float* csrValA, const int* csrRowPtrA, const int* csrColIndA,
   const float* B, int ldb, const float* beta, float* C, int ldc)
{faust_cusparseScsrmm2(handle,transA,transB,m,n,k,nnz,alpha,descrA,csrValA,csrRowPtrA,csrColIndA,B,ldb,beta,C,ldc);}

template<>
void faust_cu_csrmm<double>(cusparseHandle_t handle, 
   cusparseOperation_t transA, cusparseOperation_t transB,    
   int m, int n, int k, int nnz,       
   const double* alpha, const cusparseMatDescr_t descrA,
   const double* csrValA, const int* csrRowPtrA, const int* csrColIndA,
   const double* B, int ldb, const double* beta, double* C, int ldc)
{faust_cusparseDcsrmm2(handle,transA,transB,m,n,k,nnz,alpha,descrA,csrValA,csrRowPtrA,csrColIndA,B,ldb,beta,C,ldc);}


template<> void faust_cu_nnz<float>(cusparseHandle_t handle, cusparseDirection_t dirA, 
   int m, int n, const cusparseMatDescr_t descrA, const float* A, 
   int lda, int* nnzPerRowColumn, int* nnzTotalDevHostPtr)
{faust_cusparseSnnz(handle, dirA, m, n, descrA, A, lda, nnzPerRowColumn, nnzTotalDevHostPtr);}

template<> void faust_cu_nnz<double>(cusparseHandle_t handle, cusparseDirection_t dirA, 
   int m, int n, const cusparseMatDescr_t descrA, const double* A, 
   int lda, int* nnzPerRowColumn, int* nnzTotalDevHostPtr)
{faust_cusparseDnnz(handle, dirA, m, n, descrA, A, lda, nnzPerRowColumn, nnzTotalDevHostPtr);}


template<>
void faust_cu_dense2csr<float>(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const float* A, int lda, const int* nnzPerRow,
   float* csrValA, int* csrRowPtrA, int* csrColIndA)
{faust_cusparseSdense2csr(handle, m, n, descrA, A, lda, nnzPerRow, csrValA, csrRowPtrA, csrColIndA);}

template<>
void faust_cu_dense2csr<double>(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const double* A, int lda, const int* nnzPerRow,
   double* csrValA, int* csrRowPtrA, int* csrColIndA)
{faust_cusparseDdense2csr(handle, m, n, descrA, A, lda, nnzPerRow, csrValA, csrRowPtrA, csrColIndA);}


template<>
void faust_cu_csr2dense<float>(cusparseHandle_t handle,
 int m, int n, const cusparseMatDescr_t descrA,  
 const float* csrValA, const int* csrRowPtrA, 
 const int* csrColIndA, float* A, int lda)
{faust_cusparseScsr2dense(handle, m, n, descrA, csrValA, csrRowPtrA, csrColIndA, A, lda);}

template<>
void faust_cu_csr2dense<double>(cusparseHandle_t handle,
 int m, int n, const cusparseMatDescr_t descrA,  
 const double* csrValA, const int* csrRowPtrA, 
 const int* csrColIndA, double* A, int lda)
{faust_cusparseDcsr2dense(handle, m, n, descrA, csrValA, csrRowPtrA, csrColIndA, A, lda);}


template<>
void faust_cu_csc2dense<float>(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const float* cscValA, const int* cscRowIndA, 
   const int* cscColPtrA, float* A, int lda)
{faust_cusparseScsc2dense(handle, m, n, descrA, cscValA, cscRowIndA, cscColPtrA, A, lda);}

template<>
void faust_cu_csc2dense<double>(cusparseHandle_t handle,
   int m, int n, const cusparseMatDescr_t descrA,
   const double* cscValA, const int* cscRowIndA, 
   const int* cscColPtrA, double* A, int lda)
{faust_cusparseDcsc2dense(handle, m, n, descrA, cscValA, cscRowIndA, cscColPtrA, A, lda);}


template<>
void faust_cu_csr2csc<float>(cusparseHandle_t handle, int m, int n, int nnz,
   const float* csrVal, const int* csrRowPtr, const int* csrColInd, 
   float* cscVal, int* cscRowInd, int* cscColPtr,
   cusparseAction_t copyValues, cusparseIndexBase_t idxBase)
{faust_cusparseScsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase);}

template<>
void faust_cu_csr2csc<double>(cusparseHandle_t handle, int m, int n, int nnz,
   const double* csrVal, const int* csrRowPtr, const int* csrColInd, 
   double* cscVal, int* cscRowInd, int* cscColPtr,
   cusparseAction_t copyValues, cusparseIndexBase_t idxBase)
{faust_cusparseDcsr2csc(handle, m, n, nnz, csrVal, csrRowPtr, csrColInd, cscVal, cscRowInd, cscColPtr, copyValues, idxBase);}



