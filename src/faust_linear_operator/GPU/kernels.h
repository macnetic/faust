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

#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <unistd.h>


// for k=0 to (length-1), d_cu1[k] += d_cu2[k]
template<typename T> void kernel_add(T* d_cu1, const T* d_cu2, int length);
// for k=0 to (length-1), d_cu1[k] -= d_cu2[k]
template<typename T> void kernel_sub(T* d_cu1, const T* d_cu2, int length);
// for k=0 to (length-1), d_cu1[k] *= d_cu2[k]
template<typename T> void kernel_mult(T* d_cu1, const T* d_cu2, int length);
// for k=0 to (length-1), d_cu1[k] /= d_cu2[k]
template<typename T> void kernel_div(T* d_cu1, const T* d_cu2, int length);

// for k=0 to (length-1), d_cu1[k] += valeur
template<typename T> void kernel_add_const(T* d_cu1, T valeur, int length);
// for k=0 to (length-1), d_cu1[k] -= valeur
template<typename T> void kernel_sub_const(T* d_cu1, T valeur, int length);
// for k=0 to (length-1), d_cu1[k] *= valeur
template<typename T> void kernel_mult_const(T* d_cu1, T valeur, int length);
// for k=0 to (length-1), d_cu1[k] /= valeur
template<typename T> void kernel_div_const(T* d_cu1, T valeur, int length);

// for k=0 to (length-1), d_cu1[k] *= d_cu1[k]
template<typename T> void kernel_square(T* d_cu1, int length);
// for k=0 to (length-1), d_cu_dst[k] = d_cu1[k]*d_cu1[k]
template<typename T> void kernel_square(T* d_cu_dst, const T* d_cu_src, int length);
// for k=0 to (length-1), d_cu1[k] = sqrt(d_cu1[k])
template<typename T> void kernel_sqrt(T* d_cu1, int length);
// for k=0 to (length-1), d_cu1[k] = 1/d_cu1[k]
template<typename T> void kernel_inv(T* d_cu1, int length);
// for k=0 to (length-1), d_cu1[k] = abs(d_cu1[k]);
template<typename T> void kernel_abs(T* d_cu1, const int length);

// for k=0 to (length-1), d_cu_dst[k] = d_cu_src[k]
template<typename T, typename U> void kernel_memcpy(T* d_cu_dst, const U* d_cu_src, int length);
// for k=0 to (length-1), d_cu_dst[k] = valeur
template<typename T> void kernel_memset(T* d_cu_dst, T valeur, int length);

// convert sparse matrix (nnz, dev_src_rowind, dev_src_colind, dev_src_values) in 1-based indexing to full matrix dev_dst
// dev_dst=0; for k=0 to (nnz-1), dev_dst[(dev_src_colind[k]-1)*src_dim1+(dev_src_rowind[k]-1)] = dev_src_values[k]
template<typename T> void kernel_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1, const int src_dim2);

// add sparse matrix (nnz, dev_src_rowind, dev_src_colind, dev_src_values) in 1-based indexing to full matrix dev_dst
template<typename T> void kernel_add_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1);

// substract sparse matrix (nnz, dev_src_rowind, dev_src_colind, dev_src_values) in 1-based indexing to full matrix dev_dst
template<typename T> void kernel_sub_sparse2full(T *dev_dst, const int *dev_src_rowind, const int *dev_src_colind, const T* dev_src_values, const int nnz, const int src_dim1);

// for k=0 to (dim1-1), dst_diag[k] = src_M[k*(dim1+1)]
template<typename T> void kernel_get_diag(T* dst_diag, const T* src_M, int dim1);

// for k=0 to (dim1-1), d_cu1[k*(dim1+1)] += val
template<typename T> void kernel_add_diag_const(T* d_cu1, T val, int dim1);

// for k=0 to (dim1-1), d_cu_dst[k*(dim1+1)] = d_cu_src[k*(dim1+1)]
template<typename T> void kernel_copy_diag(T* d_cu_dst, T* d_cu_src, int dim1);

template<typename T> void kernel_set_submatrix(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_col);

// for k=0 to (length-1 )data_dst[i] = abs((data_src_th[i]-data_src_mes[i])/data_src_th[i]);
template<typename T> void kernel_relative_error(T* data_dst, const T* data_src_th, const T* data_src_mes, const int length);



#endif
