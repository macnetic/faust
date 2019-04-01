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
#ifndef __FAUST_PROX_CU_HPP__
#define __FAUST_PROX_CU_HPP__

#include <vector>
#include <iostream>
#include "algorithm"
#include "faust_gpu2cpu.h"
#include "faust_MatDense.h"
#include "faust_prox.h"

// const char * interface_prox_name="prox : ";


template<typename FPP>
void Faust::prox_sp(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_sp(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_spcol(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_spcol(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_splin(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_splin(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_normcol(Faust::MatDense<FPP,Gpu>& cu_M, FPP s)
{
    Faust::MatDense<FPP,Cpu> M;
    faust_gpu2cpu(M, cu_M);
    Faust::prox_normcol(M,s);
    cu_M = M;
}

template<typename FPP>
void Faust::prox_normlin(Faust::MatDense<FPP,Gpu> & cu_M, FPP s)
{
	Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_normlin(M,s);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_sp_pos(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_sp_pos(M,k);
   cu_M = M;
}


// M needs to be square and k must divide dimension of M
template<typename FPP>
void Faust::prox_blkdiag(Faust::MatDense<FPP,Gpu>& cu_M, int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_blkdiag(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_supp(Faust::MatDense<FPP,Gpu>& cu_M, const Faust::MatDense<FPP,Gpu>& cu_supp)
{
   Faust::MatDense<FPP,Cpu> M,supp;
   faust_gpu2cpu(M, cu_M);
   faust_gpu2cpu(supp, cu_supp);
   Faust::prox_supp(M,supp);
   cu_M = M;
}



template<typename FPP>
void Faust::prox_splincol(Faust::MatDense<FPP,Gpu>& cu_M, const faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M,supp;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_splincol(M,k);
   cu_M = M;
}
#endif


