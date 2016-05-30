#ifndef __FAUST_STOPPING_CRITERION__
#define __FAUST_STOPPING_CRITERION__

#include "faust_constant.h"
#include <iostream>


namespace Faust
{

    template<typename T>
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

          ~StoppingCriterion(){}

          bool do_continue(int current_ite, T error=-2.0)const;
          int get_crit() const{return nb_it;}
       private:
          void check_validity()const;

       private:
          // if isCriterionError then criterion is error else criterion is number of iteration
          bool  isCriterionError;
          int nb_it;   // number of iterations if !isCriterionError
          T errorThreshold;
          int maxIteration;
          // only used as stopping criterion, if isCriterionError, when error is still greater than
          static const char * class_name;
    };

}

#include "faust_StoppingCriterion.hpp"

#endif
