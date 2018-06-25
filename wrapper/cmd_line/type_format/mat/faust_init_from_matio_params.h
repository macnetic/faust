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

#ifndef __FAUST_INIT_FROM_MATIO_PARAMS_H__
#define __FAUST_INIT_FROM_MATIO_PARAMS_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

#include "faust_Params.h"
#include "faust_ParamsPalm.h"



template<typename FPP,Device DEVICE> class ParamsPalm;
template<typename FPP,Device DEVICE> class Params;
template<typename FPP,Device DEVICE> class ConstraintGeneric;

template<typename FPP,Device DEVICE>
void init_params_palm_from_matiofile(Faust::ParamsPalm<FPP,DEVICE>& params, const char* fileName, const char* variableName);

/** \brief load data matrix from ".mat file"
 * \param params
 * \param fileName
 * \param variableName
 */
template<typename FPP,Device DEVICE>
void init_params_from_matiofile(Faust::Params<FPP,DEVICE>& params, const char* fileName, const char* variableName);

template<typename FPP,Device DEVICE>
void add_constraint(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & consS,matvar_t* cons_var);

template<typename FPP,Device DEVICE>
void Display_params(Faust::Params<FPP,DEVICE> & params);

#include "faust_init_from_matio_params.hpp"
#endif
