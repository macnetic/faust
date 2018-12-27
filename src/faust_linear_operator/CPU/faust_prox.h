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
#ifndef FAUST_PROX_H
#define FAUST_PROX_H

#include "faust_MatDense.h"
#include "faust_constant.h"
#include "faust_exception.h"
#include "faust_Vect.h"

/** \brief faust_prox.h contains the projection operator: <br>
*   PALM relies on projections onto the constraint sets for each factor at each iteration, <br>
*   so the projection operator should be simple and easy to compute.
*/
namespace Faust {

    template<typename FPP>
    bool partial_sort_comp (const std::pair<int, FPP>& pair1, const std::pair<int, FPP>& pair2);

    template<typename FPP>
    void sort_idx(const std::vector<FPP> &v, std::vector<int>& idx, int s);

    template<typename FPP>
    void prox_sp(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_sp_pos(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_spcol(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splin(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_spcol_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splin_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splincol(Faust::MatDense<FPP,Cpu> &M,faust_unsigned_int k);
    template<typename FPP, typename FPP2>
    void prox_normcol(Faust::MatDense<FPP,Cpu> & M,FPP2 s);
    template<typename FPP, typename FPP2>
    void prox_normlin(Faust::MatDense<FPP,Cpu> & M,FPP2 s);
    template<typename FPP>
    void prox_supp(Faust::MatDense<FPP,Cpu> & M, const Faust::MatDense<FPP,Cpu> & supp);
//    template<typename FPP>
//    void prox_blkdiag(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);

}

#include "faust_prox.hpp"

#endif
