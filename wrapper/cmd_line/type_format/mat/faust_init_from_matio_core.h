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


#ifndef __FAUST_INIT_FROM_MATIO_CORE_H__
#define __FAUST_INIT_FROM_MATIO_CORE_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

// class Faust::Transform;
// class Faust::MatDense;
template<typename FPP,Device DEVICE> class Transform;
template<typename FPP,Device DEVICE> class MatDense;


template<typename FPP,Device DEVICE>
void init_faust_core_from_matiofile(Faust::Transform<FPP,DEVICE>& core, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void init_faust_core_from_matvar(Faust::Transform<FPP,DEVICE>& core, matvar_t* cell_var );
template<typename FPP,Device DEVICE>
void init_faust_data_from_matiofile(std::vector<Faust::MatDense<FPP,DEVICE> >& full_mat, std::vector<Faust::Transform<FPP,DEVICE> >& core, const char* fileName, const char* variableName);

#ifdef COMPILE_GPU
template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Gpu> core, const char* fileName, const char* variableName);
#endif


//template<typename FPP,Device DEVICE>
//void write_faust_spmat_list_into_matfile(const faust_spmat<FPP,DEVICE>& M, const char* fileName, const char* variableName);

template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Cpu> core, const char* fileName, const char* variableName);

template<typename FPP>
void write_faust_core_into_matvar(const Faust::Transform<FPP,Cpu> core, matvar_t** matvar, const char* variableName);





#include "faust_init_from_matio_core.hpp"

#endif
