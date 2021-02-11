/****************************************************************************/
/*                              Description:                                */
/*    file where the C++ class Faust::TransformHelper is interfaced  		*/
/*    with Matlab    														*/
/*                                                                          */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         		*/
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2021):  	Hakim HADJ-DJILANI                                  */
/*  					Nicolas Bellot, Adrien Leman, Thomas Gautrais,     	*/
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
/*		Hakim hadj-djilani : hakim.hadj-djilani@inria.fr					*/
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
#ifndef __MEX_FAUST_prox__
#define __MEX_FAUST_prox__
template <typename SCALAR, FDevice DEV>
void faust_prox(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs);
template <typename SCALAR, FDevice DEV>
void faust_prox_blockdiag(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs);
#include "faust_prox.hpp"
#endif
