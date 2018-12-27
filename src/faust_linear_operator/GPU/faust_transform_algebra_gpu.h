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
#ifndef __FAUST_Transform_ALGEBRA_GPU_H__
#define __FAUST_Transform_ALGEBRA_GPU_H__

#include "faust_constant_gpu.h"

template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class MatSparse;
template<typename FPP,Device DEVICE> class Vect;
template<typename FPP,Device DEVICE> class Transform;

// power iteration was not tasted with Faust::Transform<FPP,Gpu> because it seems that it is not tested
//template<typename FPP>
//FPP power_iteration(const Faust::Transform<FPP,Gpu> & A, const int nbr_iter_max, const FPP threshold, int & flag, cublasHandle_t cublasHandle);



//currently operator * (multiplication) can not be overloaded because GPU multiplication implies taking extra arguments cublashandle (dense matrix)  and cusparsehandle (sparse matrix)
//whereas operator must take either one or two arguments

//template<typename T>
//faust_cu_vec<T> operator*(const Faust::Transform_cu<T>& f, const faust_cu_vec<T> & v);

//template<typename T>
//faust_cu_mat<T> operator*(const Faust::Transform_cu<T>& f, const faust_cu_mat<T> & M);

#include "faust_transform_algebra_gpu.hpp"

#endif
