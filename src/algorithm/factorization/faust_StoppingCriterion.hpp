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
#ifndef __STOPPING_CRITERION_HPP__
#define __STOPPING_CRITERION_HPP__

//#include "faust_StoppingCriterion.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"

template<typename T>
const char * Faust::StoppingCriterion<T>::m_className="Faust::StoppingCriterion::";

template<typename T>
Faust::StoppingCriterion<T>::StoppingCriterion(bool isCriterionError_) : isCriterionError(isCriterionError_)
{
   if (isCriterionError_)
   {
      errorThreshold = 0.3;
      maxIteration = 10000;
   }
   else
      nb_it = 500;
}

template<typename T>
Faust::StoppingCriterion<T>::StoppingCriterion(int nb_it, bool isCriterionError, T errorThreshold, int maxIteration /* default 10000 */):
	nb_it(nb_it), isCriterionError(isCriterionError), errorThreshold(errorThreshold), maxIteration(maxIteration)
{
	check_validity();
}


template<typename T>
void Faust::StoppingCriterion<T>::check_validity()const
{
   if (isCriterionError)
   {
      if (errorThreshold>1 || maxIteration < 0)
      {
        handleError(m_className,"check_validity : errorThreshold must be strictly greater than 1 and maxIteration must be strictly positive");
      }
   }
   else if (nb_it < 0)
   {
     handleError(m_className,"::check_validity : nb_it must be positive");
   }
}

// current_ite in zero-based indexing
template<typename T>
bool Faust::StoppingCriterion<T>::do_continue(int current_ite, T current_error /* = -2.0 */)const
{

   if (!isCriterionError) // if criterion is number of iteration, current_error does not matter
      return current_ite<nb_it ? true : false;
   else if (isCriterionError && current_error != -2.0)
      if (current_error < errorThreshold)
         return false;
      else if (current_ite <  maxIteration) // and current_error >= errorThreshold
         return true;
      else // if current_error >= errorThreshold and current_ite >= maxIteration
      {
         std::cerr << "warning in Faust::StoppingCriterion<T>::do_continue : number of maximum iterations has been reached and current error is still greater than the threshold" << std::endl;
         return true;
      }
   else // if criterion is error and current_error has not been initialized
   {
     handleError(m_className,"check_validity : when stopping criterion is error, the current error needs to be given as second parameter");
   }
}

template<typename T>
void Faust::StoppingCriterion<T>::Display() const
{
	std::cout << "StoppingCriterion obj:" << std::endl;
	std::cout << "\tnb_it="<< nb_it << std::endl;
	std::cout << "\tisCriterionError="<< isCriterionError << std::endl;
	std::cout << "\terrorTreshold="<< errorThreshold << std::endl;
	std::cout << "\tmaxIteration="<< maxIteration << std::endl;
}
#endif
