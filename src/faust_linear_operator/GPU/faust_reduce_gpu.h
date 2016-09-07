/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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

#ifndef __FAUST_CU_REDUCE_H__
#define __FAUST_CU_REDUCE_H__

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>

/*template<typename faust_real> class faust_cu_vec;
template<typename faust_real> class faust_cu_matrix;

template<typename faust_real>
faust_real faust_cu_reduce(const faust_cu_vec<faust_real>& v);

template<typename faust_real>
faust_real faust_cu_reduce(const faust_cu_matrix<faust_real>& M);*/


template<typename FPP>
FPP faust_cu_sum(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_max(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_min(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_norm(const FPP* data, const int nb_el);



#endif
