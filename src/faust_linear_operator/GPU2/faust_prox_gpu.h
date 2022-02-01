/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __PROX_GPU_H__
#define __PROX_GPU_H__

#include <vector>
#include <iostream>
#include <algorithm>
#include "faust_prox.h"
#include "faust_constant.h"
#include "faust_MatDense_gpu.h"
#include "faust_exception.h"

namespace Faust
{
	template<typename FPP>
		void prox_sp(MatDense<FPP,GPU2> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false, const bool pure_gpu=true);
	template<typename FPP>
		void prox_sp_pos(MatDense<FPP, GPU2> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false );
	template<typename FPP>
		void prox_spcol(MatDense<FPP,GPU2> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false, const bool pure_gpu=true);
	template<typename FPP>
		void prox_splin(MatDense<FPP,GPU2> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false, const bool pure_gpu=true);
	template<typename FPP>
		void prox_splincol(MatDense<FPP,GPU2> &M,faust_unsigned_int k, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_supp(MatDense<FPP,GPU2> & M, const MatDense<FPP,GPU2> & supp, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_const(MatDense<FPP,GPU2> & M, const MatDense<FPP,GPU2> & supp, const bool normalized=false, const bool pos=false);
	template<typename FPP>
		void prox_id(MatDense<FPP,GPU2> & M, const bool normalized=false, const bool pos=false);

	template<typename FPP, typename FPP2>
		void prox_normcol(MatDense<FPP,GPU2> & M,FPP2 s, const bool normalized=false, const bool pos=false);
	template<typename FPP, typename FPP2>
		void prox_normlin(MatDense<FPP,GPU2> & M,FPP2 s, const bool normalized=false, const bool pos=false);
	//    template<typename FPP>
	//    void prox_blkdiag(MatDense<FPP,GPU2> & M,faust_unsigned_int k);
	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,GPU2> & M, std::vector<faust_unsigned_int>& m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized = false, const bool pos = false);
	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,GPU2> & M, MatDense<FPP,GPU2> mn_vec, const bool normalized=true, const bool pos=false);

	template<typename FPP> void pre_prox_pos(MatDense<FPP,GPU2> & M);
	template<typename FPP> void prox_hankel(MatDense<FPP, GPU2> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP> void prox_toeplitz(MatDense<FPP, GPU2> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP> void prox_circ(MatDense<FPP, GPU2> & M, const bool normalized = true, const bool pos = false);

	template<typename FPP> void prox_skperm(MatDense<FPP, GPU2> & M,const unsigned int k,  const bool normalized=true, const bool pos=false);

}
#include "faust_prox_gpu.hpp"

#endif
