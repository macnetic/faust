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
#ifndef __FAUST_Transform_ALGEBRA_HPP__
#define __FAUST_Transform_ALGEBRA_HPP__

#include <iostream>
#include "faust_linear_algebra.h"
#include "faust_MatDense.h"
#include "faust_Vect.h"
#include <Eigen/QR>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

#include "faust_MatSparse.h"
#include "faust_Transform.h"
#include "faust_exception.h"
#include <complex>

	//////////FONCTION Faust::MatDense<FPP,Cpu> - Faust::MatDense<FPP,Cpu> ////////////////////

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif




// const char * core_algebra_name="Faust::Transform<FPP,Cpu>_algebra : ";










//////////// modif AL AL
 
//////////// modif AL AL

template<typename FPP>
#ifdef __COMPILE_TIMERS__
	Faust::Vect<FPP,Cpu> Faust::operator*(Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v)
#else		
	Faust::Vect<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v)
#endif
{

	return f.multiply(v);
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M)
{
	return f.multiply(M);
}

#endif
