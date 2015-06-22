#ifndef __FAUST_STOPPING_CRITERION__
#define __FAUST_STOPPING_CRITERION__

#include "faust_constant.h"
#include <iostream>

class stopping_criterion
{
   public:
      stopping_criterion():
         isCriterionError(false),
         nb_it(500){}

      stopping_criterion(int nb_it_):
         isCriterionError(false),nb_it(nb_it_){check_validity();}
         
      stopping_criterion(faust_real errorThreshold_, int maxIteration_=10000):
         isCriterionError(true),errorThreshold(errorThreshold_),maxIteration(maxIteration_){check_validity();}

      stopping_criterion(bool isCriterionError_);

      ~stopping_criterion(){}

      bool do_continue(int current_ite, faust_real error=-2.0)const;

   private:
      void check_validity()const;

   private:
      int nb_it;   // number of iterations if !isCriterionError
      // if isCriterionError then criterion is error else criterion is number of iteration
      bool  isCriterionError;
      faust_real errorThreshold;
      // only used as stopping criterion, if isCriterionError, when error is still greater than 
      int maxIteration;
};


#endif
