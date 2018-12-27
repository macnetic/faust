/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2018):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#include "faust_constant.h"

#ifdef __GEMM_WITH_OPENBLAS__
#include "faust_cblas_algebra.h"
#include <complex>
#include "faust_exception.h"

namespace Faust
{

    ////  void Faust::cblas_gemm(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const T alpha, const T* Adata,const int lda,const T* Bdata,const int ldb, const T beta, T* Cdata, const int ldc);
    template<> void cblas_gemm<float>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const float alpha, const float* Adata,const int lda,const float* Bdata,const int ldb, const float beta, float* Cdata, const int ldc)
    {
        cblas_sgemm(order,transa,transb,dim1,dim2,nbColOpA,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc);
        #ifdef	FAUST_VERBOSE
            std::cout<<"*** inside float Faust::cblas_gemm"<<std::endl;
        #endif
    }

    template<> void cblas_gemm<double>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const double alpha, const double* Adata,const int lda,const double* Bdata,const int ldb, const double beta, double* Cdata, const int ldc)
    {
        cblas_dgemm(order,transa,transb,dim1,dim2,nbColOpA,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc);
        #ifdef	FAUST_VERBOSE
            std::cout<<"*** inside double Faust::cblas_gemm"<<std::endl;
        #endif

    }
	
    
    template<> void cblas_gemm<std::complex<float> >(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const std::complex<float> alpha, const std::complex<float>* Adata,const int lda,const std::complex<float>* Bdata,const int ldb, const std::complex<float> beta, std::complex<float>* Cdata, const int ldc)
    {
        handleError("cblas_gemm<std::complex<float>>","not yet implemented");
    }


    template<> void cblas_gemm<std::complex<double> >(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const std::complex<double> alpha, const std::complex<double>* Adata,const int lda,const std::complex<double>* Bdata,const int ldb, const std::complex<double> beta, std::complex<double>* Cdata, const int ldc)
    {
        handleError("cblas_gemm<std::complex<double>>","not yet implemented");
    }		



    template<> void cblas_gemv<float>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const float alpha,const float* Adata,const int lda,const float* Xdata,const int incX,const float beta,float* Ydata,const int incY)
    {
        cblas_sgemv(order,transa,dim1,dim2,alpha,Adata,lda,Xdata,incX,beta,Ydata,incY);
        #ifdef	FAUST_VERBOSE
            std::cout<<"***inside float Faust::cblas_gemv"<<std::endl;
        #endif
    }

    template<> void cblas_gemv<double>(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const double alpha,const double* Adata,const int lda,const double* Xdata,const int incX,const double beta,double* Ydata,const int incY)
    {
        cblas_dgemv(order,transa,dim1,dim2,alpha,Adata,lda,Xdata,incX,beta,Ydata,incY);
        #ifdef	FAUST_VERBOSE
            std::cout<<"***inside double Faust::cblas_gemv"<<std::endl;
        #endif
    }


       template<> void cblas_gemv<std::complex<double> >(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const std::complex<double> alpha,const std::complex<double>* Adata,const int lda,const std::complex<double>* Xdata,const int incX,const std::complex<double> beta,std::complex<double>* Ydata,const int incY)
    {
        handleError("cblas_gemv<std::complex<double>>","not yet implemented");
    }

           template<> void cblas_gemv<std::complex<float> >(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const std::complex<float> alpha,const std::complex<float>* Adata,const int lda,const std::complex<float>* Xdata,const int incX,const std::complex<float> beta,std::complex<float>* Ydata,const int incY)
    {
        handleError("cblas_gemv<std::complex<float>>","not yet implemented");
    }		

}

#endif
