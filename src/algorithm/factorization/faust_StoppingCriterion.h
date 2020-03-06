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
#ifndef __FAUST_STOPPING_CRITERION__
#define __FAUST_STOPPING_CRITERION__

#include "faust_constant.h"
#include <iostream>


namespace Faust
{

    template<typename T = double>
    class StoppingCriterion
    {
       public:
          StoppingCriterion():
             isCriterionError(false),
             nb_it(500){}

          StoppingCriterion(int nb_it_):
             isCriterionError(false),nb_it(nb_it_){check_validity();}

          StoppingCriterion(T errorThreshold_, int maxIteration_=10000):
             isCriterionError(true),errorThreshold(errorThreshold_),maxIteration(maxIteration_){check_validity();}

          StoppingCriterion(bool isCriterionError_);

          StoppingCriterion(int nb_it, bool isCriterionError, T errorThreshold, int maxIteration=10000);


          ~StoppingCriterion(){}

          bool do_continue(int current_ite, T error=-2.0)const;
          int get_crit() const{return nb_it;}
		  bool isCriterionErr() const {return isCriterionError;}
		  void Display() const;
		  static const T NO_ERROR_PASSED;
       private:
          void check_validity()const;

       private:
          // if isCriterionError then criterion is error else criterion is number of iteration
          bool  isCriterionError;
          int nb_it;   // number of iterations if !isCriterionError
          T errorThreshold;
          int maxIteration;
          // only used as stopping criterion, if isCriterionError, when error is still greater than
          static const char * m_className;
    };

}

#include "faust_StoppingCriterion.hpp"

#endif
