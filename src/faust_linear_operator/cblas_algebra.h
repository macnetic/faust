#ifdef __GEMM_WITH_OPENBLAS__
#ifndef CBLAS_ALGEBRA_H
#define CBLAS_ALGEBRA_H

#include "cblas.h"
#include <iostream>

/*! \brief cblas_gemm call the function cblas_sgemm or cblas_dgemm, depend to the precision (double or float), from the external library openBLAS    <br> 
* This function is used to multiply the matrices.
*/
template<typename T>
void cblas_gemm(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const T alpha, const T* Adata,const int lda,const T* Bdata,const int ldb, const T beta, T* Cdata, const int ldc);

// cblas_sgemv(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);


/*! \brief cblas_gemv call the function cblas_sgemv or cblas_dgemv, depend to the precision (double or float), from the external library openBLAS    <br> 
* This function is used to multiply the matrices.
*/
template<typename T>
void cblas_gemv(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const T alpha,const T* Adata,const int lda,const T* Xdata,const int incX,const T beta,T* Ydata,const int incY);

#endif
#endif


