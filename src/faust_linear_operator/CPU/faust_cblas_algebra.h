/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifdef __GEMM_WITH_OPENBLAS__
#ifndef CBLAS_ALGEBRA_H
#define CBLAS_ALGEBRA_H

#include "cblas.h"
#include <iostream>

namespace Faust
{

    /*! \fn Faust::cblas_gemm
    *   \brief call the function cblas_sgemm or cblas_dgemm, depend to the precision (double or float), from the external library openBLAS    <br>
    * This function is used to multiply the matrices.
    */
    template<typename T>
    void cblas_gemm(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int dim1,const int dim2,const int nbColOpA, const T alpha, const T* Adata,const int lda,const T* Bdata,const int ldb, const T beta, T* Cdata, const int ldc);

// cblas_sgemv(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);


    /*! \fn Faust::cblas_gemv
    *   \brief call the function cblas_sgemv or cblas_dgemv, depend to the precision (double or float), from the external library openBLAS    <br>
    * This function is used to multiply the matrices.
    */
    template<typename T>
    void cblas_gemv(const CBLAS_ORDER order,const CBLAS_TRANSPOSE transa,const int dim1,const int dim2,const T alpha,const T* Adata,const int lda,const T* Xdata,const int incX,const T beta,T* Ydata,const int incY);
}

#endif
#endif


