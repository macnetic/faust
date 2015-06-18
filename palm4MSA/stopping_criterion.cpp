#include "stopping_criterion.h"
#include <iostream>
#include <cstdlib>

stopping_criterion::stopping_criterion(bool isCriterionError_) : isCriterionError(isCriterionError_)
{
   if (isCriterionError_)
   {
      errorThreshold = 0.3;
      maxIteration = 10000;
   }
   else
      nb_it = 500;
}

void stopping_criterion::check_validity()const
{
   if (isCriterionError)
      if (errorThreshold>1 || maxIteration < 0)
      {
         std::cerr << "error in stopping_criterion::check_validity" << std::endl;
         exit(EXIT_FAILURE);
      }
   else if (nb_it < 0) 
   {
      std::cerr << "error in stopping_criterion::check_validity" << std::endl;
      exit(EXIT_FAILURE);
   }
}

// current_ite in zero-based indexing
bool stopping_criterion::do_continue(int current_ite, faust_real current_error /* = -2.0 */)const
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
         std::cerr << "warning in stopping_criterion::do_continue : number of maximum iterations has been reached and current error is still greater than the threshold" << std::endl;
         return true;
      }
   else // if criterion is error and current_error has not been initialized
   {
      std::cerr << "error in stopping_criterion::check_validity : when stopping criterion is error, the current error needs to be given as second parameter" << std::endl;
      exit(EXIT_FAILURE);
   }
}

