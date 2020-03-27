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
#ifndef FAUST_Transform_ALGEBRA_H
#define FAUST_Transform_ALGEBRA_H

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_linear_algebra.h"
#include <complex>


namespace Faust
{

	template<typename FPP,FDevice DEVICE> class MatDense;
	template<typename FPP,FDevice DEVICE> class Vect;
	template<typename FPP,FDevice DEVICE> class Transform;


	// compute absolute value of complex scalar
	template<typename FPP>
	FPP absValue(std::complex<FPP> cplxScalar){return std::abs(cplxScalar);}
	
	// compute absolute value of real scalar
	template<typename FPP>
	FPP absValue(FPP realScalar){return fabs(realScalar);}
	
	// if measure of time is down Faust::Transform<FPP,Cpu> is no longer constant during multiplication because a measure of time is an attribute to the faust::Transform
	template<typename FPP>
	#ifdef __COMPILE_TIMERS__
		Faust::Vect<FPP,Cpu> operator*(Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v);
	#else		
		Faust::Vect<FPP,Cpu> operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v);
	#endif

	template<typename FPP>
	Faust::MatDense<FPP,Cpu> operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu> & M);

}

#include "faust_transform_algebra.hpp"

#endif
