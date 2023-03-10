/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2021):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef FAUST_PROX_GEN_H
#define FAUST_PROX_GEN_H

#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_constant.h"
#include "faust_Vect.h"
#include "faust_prox.h"
#ifdef USE_GPU_MOD
#include "faust_prox_gpu.h"
#endif
#include <iostream>

/** \brief faust_prox_gen.h defines the projection operators exactly as faust_prox.h except that this module is doubly transparent:
 * 1) it returns MatGeneric-s (type concrete type, MatSparse or MatDense, is automatically decided in order to optimize the memory footprint).
 * 2) It works for FDevice GPU2 and Cpu.<br>
*/
namespace Faust
{


	template<typename FPP, FDevice DEV> faust_unsigned_int sparse_size(faust_unsigned_int nnz, faust_unsigned_int nrows);
	template<typename FPP, FDevice DEV> faust_unsigned_int dense_size(faust_unsigned_int nrows, faust_unsigned_int ncols);

    /**
     * Decides which output format to use when appliying the SP prox. op. to M. Either M or spM is the output, it depends on the byte sizes. The minimum memory fingerprint is targeted.
     *
     *  \param M: the input matrix to project (which could be the output too, if dense format is chosen).
     *  \param k: the sparsity parameter (pseudo-)norm_0 of the output matrix.
     *  \param normalized: true to normalize the output matrix.
     *  \param pos: true to filter negative values of M before applying the de prox.
     *  \param forcedType: used to choose explicitely the output format with values Sparse or Dense (MatSparse or MatDense).
     * \return the prox image as a MatGeneric matrix (which can be a MatSparse or a MatDense). if the returned matrix is a MatDense, the memory spaced used is the same as the input (that is M), otherwise the matrix is a MatSparse and is allocated in the heap memory (and must be freed by the callee).
     */
	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_sp_gen(MatDense<FPP,DEV> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_skperm_gen(MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized=true, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_splincol_gen(MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized=true, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_spcol_gen(MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized=true, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_splin_gen(MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized=true, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV, typename FPP2>
		MatGeneric<FPP,DEV>* prox_normcol_gen(MatDense<FPP,DEV> & M,FPP2 s, const bool normalized=false, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV, typename FPP2>
		MatGeneric<FPP,DEV>* prox_normlin_gen(MatDense<FPP,DEV> & M,FPP2 s, const bool normalized=false, const bool pos=false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_blockdiag_gen(MatDense<FPP,DEV> & M, std::vector<faust_unsigned_int>& m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized = false, const bool pos = false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_blockdiag_gen(MatDense<FPP,DEV> & M, MatDense<FPP,DEV> mn_vec, const bool normalized=false, const bool pos=false, const MatType forcedType=None);
	template<typename FPP, FDevice DEV> MatGeneric<FPP,DEV>* prox_hankel_gen(MatDense<FPP, DEV> & M, const bool normalized = true, const bool pos = false, const MatType forcedType=None);
	template<typename FPP, FDevice DEV> MatGeneric<FPP,DEV>* prox_toeplitz_gen(MatDense<FPP, DEV> & M, const bool normalized = true, const bool pos = false, const MatType forcedType=None);
	template<typename FPP, FDevice DEV> MatGeneric<FPP,DEV>* prox_circ_gen(MatDense<FPP, DEV> & M, const bool normalized = true, const bool pos = false, const MatType forcedType=None);
	template<typename FPP, FDevice DEV> MatGeneric<FPP,DEV>* prox_anticirc_gen(MatDense<FPP, DEV> & M, const bool normalized = true, const bool pos = false, const MatType forcedType=None);

	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_supp_gen(MatDense<FPP,DEV> & M, const MatDense<FPP,DEV> & supp, const bool normalized=true, const bool pos=false,  const MatType forcedType=None);
	template<typename FPP, FDevice DEV>
		MatGeneric<FPP,DEV>* prox_const_gen(MatDense<FPP,DEV> & M, const MatDense<FPP,DEV> & supp, const bool normalized=true, const bool pos=false,  const MatType forcedType=None);
	template<typename FPP, FDevice DEV>
		Faust::MatGeneric<FPP,DEV>* prox_id_gen(Faust::MatDense<FPP,DEV> & M, const bool normalized=false, const bool pos=false,  const MatType forcedType=None);

  template<typename FPP, FDevice DEV>
  MatGeneric<FPP,DEV>* prox_triu_sp_gen(MatDense<FPP,DEV> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false, const MatType forcedType=None);
  template<typename FPP, FDevice DEV>
  MatGeneric<FPP,DEV>* prox_tril_sp_gen(MatDense<FPP,DEV> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false, const MatType forcedType=None);

}

#include "faust_prox_gen.hpp"

#endif
