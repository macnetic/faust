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



#ifndef __FAUST_INIT_FROM_MATIO_MAT_H__
#define __FAUST_INIT_FROM_MATIO_MAT_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include "faust_constant.h"
#include <vector>
#include "faust_MatDense.h"
#include "faust_MatSparse.h"

template<typename FPP,Device DEVICE> class MatSparse;
template<typename FPP,Device DEVICE> class MatDense;


template<typename FPP,Device DEVICE>
void init_faust_mat_from_matio(Faust::MatDense<FPP,DEVICE>& M, const char* fileName, const char* variableName);

template<typename FPP,Device DEVICE>
void init_faust_spmat_from_matio(Faust::MatSparse<FPP,DEVICE>& S, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void write_faust_mat_list_into_matfile(const std::vector< Faust::MatDense<FPP,DEVICE> >& M, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void init_faust_mat_vector_from_matiofile( std::vector<Faust::MatDense<FPP,DEVICE> > & vec_M, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void init_mat_from_matvar(Faust::MatDense<FPP,DEVICE> & M,matvar_t** var);
template<typename FPP,Device DEVICE>
void init_spmat_from_matvar(Faust::MatSparse<FPP,DEVICE> & M,matvar_t* var);






template<typename FPP,Device DEVICE>
void write_faust_mat_list_into_matvar(const std::vector<Faust::MatDense<FPP,DEVICE> >& M,matvar_t** matvar, const char* variableName);


//passer l adresse du pointeur en parametre, pas un pointeur de pointeur pour matvar_t** matvar
template<typename FPP,Device DEVICE>
void write_faust_mat_into_matvar(const Faust::MatDense<FPP,DEVICE>& M,matvar_t** matvar, const char* variableName);
template<typename FPP,Device DEVICE>
void write_faust_mat_into_matfile(const Faust::MatDense<FPP,DEVICE>& M, const char* fileName, const char* variableName);

/// A MODIFIER : pour l instant les matrices creuse sont stockés en matrices dense dans matlab
template<typename FPP,Device DEVICE>
void write_faust_spmat_list_into_matfile(const Faust::MatSparse<FPP,DEVICE>& M, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void write_faust_spmat_list_into_matvar(const std::vector<Faust::MatSparse<FPP,DEVICE> >& M,matvar_t** matvar, const char* variableName);





#include "faust_init_from_matio_mat.hpp"

#endif
