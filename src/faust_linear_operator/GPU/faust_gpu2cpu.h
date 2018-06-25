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
#ifndef __FAUST_GPU2CPU_H__
#define __FAUST_GPU2CPU_H__
#include "faust_cuda.h"


//modif AL AL
#include "faust_constant_gpu.h"
//#include "faust_MatDense_gpu.h"


// modif AL AL
template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class Vect;

#ifdef __COMPILE_SPMAT__
    template<typename FPP,Device DEVICE> class MatSparse;
    template<typename FPP,Device DEVICE> class Transform;
#endif

//! \brief faust_gpu2cpu copy an CPU vector/matrix into a GPU vector/matrix depend of the overloading function.
template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Vect<FPP,Cpu>& v, const Faust::Vect<FPP1,Gpu>& cu_v, cudaStream_t stream=0);
template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::MatDense<FPP,Cpu>& M, const Faust::MatDense<FPP1,Gpu>& cu_M, cudaStream_t stream=0);
#ifdef __COMPILE_SPMAT__
    template<typename FPP, typename FPP1>
    void faust_gpu2cpu(Faust::MatSparse<FPP,Cpu>& S, const Faust::MatSparse<FPP1,Gpu>& cu_S, cudaStream_t stream=0);
    template<typename FPP, typename FPP1>
    void faust_gpu2cpu(Faust::Transform<FPP,Cpu>& fcore, const Faust::Transform<FPP1,Gpu>& cu_fcore);
#endif


#include "faust_gpu2cpu.hpp"

#endif
