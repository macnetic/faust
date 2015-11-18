


#ifdef __GEMM_WITH_OPENBLAS__

#ifndef CBLAS_ALGEBRA_H
#define CBLAS_ALGEBRA_H

#include "cblas.h"
#include <iostream>

template<typename T>
void cblas_gemm(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const T alpha, const T* Adata,const int lda,const T* Bdata,const int ldb, const T beta, T* Cdata, const int ldc);



// cblas_sgemv(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);
template<typename T>
void cblas_gemv(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const T alpha,const T* Adata,const int lda,const T* Xdata,const int incX,const T beta,T* Ydata,const int incY);





#endif


#endif