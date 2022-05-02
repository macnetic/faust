/*
* Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above
*    copyright notice, this list of conditions and the following
*    disclaimer in the documentation and/or other materials provided
*    with the distribution.
*
* 3. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* */

// The code is from the scipy package, sligthly modified in order to use different template types
// https://github.com/scipy/scipy/blob/main/scipy/sparse/sparsetools/csr.h

#ifndef __FAUST_SCIPY__
#define __FAUST_SCIPY__

/*
 * Slice columns given as an array of indices (pass 1).
 * This pass counts idx entries and computes a new indptr.
 *
 * Input Arguments:
 *   I  n_idx           - number of indices to slice
 *   I  col_idxs[n_idx] - indices to slice
 *   I  n_row           - major axis dimension
 *   I  n_col           - minor axis dimension
 *   I  Ap[n_row+1]     - indptr
 *   I  Aj[nnz(A)]      - indices
 *
 * Output Arguments:
 *   I  col_offsets[n_col] - cumsum of index repeats
 *   I  Bp[n_row+1]        - new indptr
 *
 */
template<class I, class J>
void csr_column_index1(const J n_idx,
                       const J col_idxs[],
                       const J n_row,
                       const J n_col,
                       const I Ap[],
                       const I Aj[],
                       I col_offsets[],
                       I Bp[]);


/*
 * Slice columns given as an array of indices (pass 2).
 * This pass populates indices/data entries for selected columns.
 *
 * Input Arguments:
 *   I  col_order[n_idx]   - order of col indices
 *   I  col_offsets[n_col] - cumsum of col index counts
 *   I  nnz                - nnz(A)
 *   I  Aj[nnz(A)]         - column indices
 *   T  Ax[nnz(A)]         - data
 *
 * Output Arguments:
 *   I  Bj[nnz(B)] - new column indices
 *   T  Bx[nnz(B)] - new data
 *
 */
template<class I, class J, class T>
void csr_column_index2(const J col_order[],
                       const I col_offsets[],
                       const J nnz,
                       const I Aj[],
                       const T Ax[],
                       I Bj[],
                       T Bx[]);


#include "faust_scipy.hpp"
#endif
