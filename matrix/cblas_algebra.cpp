#include "faust_constant.h"

#ifdef __GEMM_WITH_OPENBLAS__

#include "cblas_algebra.h" 


// void cblas_gemm(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const T alpha, const T* Adata,const int lda,const T* Bdata,const int ldb, const T beta, T* Cdata, const int ldc);
template<> void cblas_gemm<float>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const float alpha, const float* Adata,const int lda,const float* Bdata,const int ldb, const float beta, float* Cdata, const int ldc)
{
	cblas_sgemm(order,transa,transb,dim1,dim2,nbColOpA,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc);
	std::cout<<"*** inside float cblas_gemm"<<std::endl;
}




template<> void cblas_gemm<double>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const double alpha, const double* Adata,const int lda,const double* Bdata,const int ldb, const double beta, double* Cdata, const int ldc)
{
	cblas_dgemm(order,transa,transb,dim1,dim2,nbColOpA,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc);
	std::cout<<"*** inside double cblas_gemm"<<std::endl;
}




template<> void cblas_gemv<float>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const float alpha,const float* Adata,const int lda,const float* Xdata,const int incX,const float beta,float* Ydata,const int incY)
{
	cblas_sgemv(order,transa,dim1,dim2,alpha,Adata,lda,Xdata,incX,beta,Ydata,incY);
	std::cout<<"***inside float cblas_gemv"<<std::endl;
}

template<> void cblas_gemv<double>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const double alpha,const double* Adata,const int lda,const double* Xdata,const int incX,const double beta,double* Ydata,const int incY)
{
	cblas_dgemv(order,transa,dim1,dim2,alpha,Adata,lda,Xdata,incX,beta,Ydata,incY);
	std::cout<<"***inside double cblas_gemv"<<std::endl;
}


#endif