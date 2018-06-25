/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
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
#ifndef __PROX_GPU_H__
#define __PROX_GPU_H__

#include "faust_MatDense_gpu.h"
#include "faust_constant_gpu.h"
#include "faust_exception.h"

namespace Faust
{

    template<typename FPP>
    void prox_sp(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_sp_pos(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_spcol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splin(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splincol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);


    template<typename FPP>
    void prox_normcol(Faust::MatDense<FPP,Gpu> & M,FPP s);
    template<typename FPP>
    void prox_normlin(Faust::MatDense<FPP,Gpu> & M,FPP s);
    template<typename FPP>
    void prox_supp(Faust::MatDense<FPP,Gpu> & M, const Faust::MatDense<FPP,Gpu> & supp);
    template<typename FPP>
    void prox_blkdiag(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
    template<typename FPP>
    void prox_splincol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);

}


#include "faust_prox_gpu.hpp"

#endif
