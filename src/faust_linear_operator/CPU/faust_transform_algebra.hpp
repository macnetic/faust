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
#ifndef __FAUST_Transform_ALGEBRA_HPP__
#define __FAUST_Transform_ALGEBRA_HPP__

#include <iostream>
#include "faust_linear_algebra.h"
#include "faust_MatDense.h"
#include "faust_Vect.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "faust_MatSparse.h"
#include "faust_Transform.h"
#include "faust_exception.h"

	//////////FONCTION Faust::MatDense<FPP,Cpu> - Faust::MatDense<FPP,Cpu> ////////////////////

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif



// const char * core_algebra_name="Faust::Transform<FPP,Cpu>_algebra : ";
template<typename FPP>
FPP Faust::power_iteration(const  Faust::Transform<FPP,Cpu> & A, const int nbr_iter_max,FPP threshold, int & flag)
{


 const int nb_col = A.getNbCol();
   int i = 0;
   flag = 0;

   if (nbr_iter_max <= 0)
      handleError("faust_transform_algebra "," power_iteration :  nbr_iter_max <= 0");
   if (nb_col != A.getNbRow())
      handleError("faust_transform_algebra "," power_iteration : Faust::Transform_gpu<FPP,Cpu> 1 must be a squared matrix");

   Faust::Vect<FPP,Cpu> xk(nb_col);
   xk.setOnes();
   Faust::Vect<FPP,Cpu> xk_norm(nb_col);
   FPP lambda_old=1.0;
   FPP lambda = 0.0;
   while(fabs(lambda_old-lambda)>threshold && i<nbr_iter_max)
   {
      i++;
      lambda_old = lambda;
      xk_norm = xk;
      xk_norm.normalize();
      xk = A*xk_norm;
      lambda = xk_norm.dot(xk);
   }
   flag = (i<nbr_iter_max)?i:-1;
   return lambda;










}

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
