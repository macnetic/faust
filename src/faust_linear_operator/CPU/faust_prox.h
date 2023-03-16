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
#ifndef FAUST_PROX_H
#define FAUST_PROX_H

#include "faust_MatDense.h"
#include "faust_constant.h"
#include "faust_exception.h"
#include "faust_Vect.h"
#include <Eigen/Core>
#include <queue>
#include <iostream>
#include <cassert>
#include <utility>
#include <limits>
#include <cstdlib>
#include <type_traits>
#include <cmath>

/** \brief faust_prox.h contains the projection operator: <br>
*   PALM relies on projections onto the constraint sets for each factor at each iteration, <br>
*   so the projection operator should be simple and easy to compute.
*/
namespace Faust
{

	template<typename FPP>
		bool partial_sort_comp (const std::pair<int, FPP>& pair1, const std::pair<int, FPP>& pair2);

	template<typename FPP>
		void sort_idx(const std::vector<FPP> &v, std::vector<int>& idx, int s);

	template<typename FPP>
		void prox_sp(MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_sp_pos(MatDense<FPP, Cpu> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false );
	template<typename FPP>
		void prox_spcol(MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_splin(MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_splincol(MatDense<FPP,Cpu> &M,faust_unsigned_int k, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_supp(MatDense<FPP,Cpu> & M, const MatDense<FPP,Cpu> & supp, const bool normalized=true, const bool pos=false);
	template<typename FPP>
		void prox_const(MatDense<FPP,Cpu> & M, const MatDense<FPP,Cpu> & supp, const bool normalized=false, const bool pos=false);

	template<typename FPP, typename FPP2>
		void prox_normcol(MatDense<FPP,Cpu> & M,FPP2 s, const bool normalized=false, const bool pos=false);
	template<typename FPP, typename FPP2>
		void prox_normlin(MatDense<FPP,Cpu> & M,FPP2 s, const bool normalized=false, const bool pos=false);
	//    template<typename FPP>
	//    void prox_blkdiag(MatDense<FPP,Cpu> & M,faust_unsigned_int k);
	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,Cpu> & M, std::vector<faust_unsigned_int>& m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized = true, const bool pos = false);
	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,Cpu> & M, MatDense<FPP,Cpu> mn_vec, const bool normalized=true, const bool pos=false);

	template<typename FPP> void pre_prox_pos(MatDense<FPP,Cpu> & M);
	template<typename FPP> void prox_hankel(MatDense<FPP, Cpu> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP> void prox_toeplitz(MatDense<FPP, Cpu> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP> void prox_circ(MatDense<FPP, Cpu> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP> void prox_anticirc(MatDense<FPP, Cpu> & M, const bool normalized = true, const bool pos = false);
	template<typename FPP>
		void prox_skperm(MatDense<FPP, Cpu> & M,const unsigned int k,  const bool normalized=true, const bool pos=false);

	template<typename FPP> faust_unsigned_int sparse_size(faust_unsigned_int nnz, faust_unsigned_int nrows);
	template<typename FPP> faust_unsigned_int dense_size(faust_unsigned_int nrows, faust_unsigned_int ncols);
	template<typename FPP>
		void prox_id(MatDense<FPP,Cpu> & M, const bool normalized=false, const bool pos=false);

  template<typename FPP>
  void prox_tri_sp(MatDense<FPP, Cpu> & M, faust_unsigned_int k, bool upper, const bool normalized=true, const bool pos=false);
  template<typename FPP>
  void prox_triu_sp(MatDense<FPP, Cpu> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false);
  template<typename FPP>
  void prox_tril_sp(MatDense<FPP, Cpu> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false);
  template<typename FPP>
  void prox_symm_sp(MatDense<FPP, Cpu> & M, faust_unsigned_int k, const bool normalized=true, const bool pos=false);

}

#include "faust_prox.hpp"

#endif
