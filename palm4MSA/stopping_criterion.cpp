#include "stopping_criterion.h"

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

void stopping_criterion::check_validity()
{
   if (isCriterionError)
      if (errorThreshold>1 || maxIteration < 0)
      {
         cerr << "error in stopping_criterion::check_validity" << endl;
         exit(EXIT_FAILURE);
      }
   else if (nb_it < 0) 
   {
      cerr << "error in stopping_criterion::check_validity" << endl;
      exit(EXIT_FAILURE);
   }
}

// current_ite in zero-based indexing
bool stopping_criterion::do_continue(int current_ite, faust_real current_error /* = -2.0 */)
{
   if (isCriterionError && current_error != -2.0) 
      if (current_error < errorThreshold)
         return false;
      else if (current_ite <  maxIteration) // and current_error >= errorThreshold
         return true;
      else // if current_error >= errorThreshold and current_ite >= maxIteration
      {
         cerr << "warning in stopping_criterion::do_continue : number of maximum iterations has been reached and current error is still greater than the threshold" <<endl;
         return true;
      }
      else // if current_error>=errorThreshold and 
         return false;

   else if (current_error == -2.0)
   {
      cerr << "error in stopping_criterion::check_validity : when stopping criterion is error, the current error needs to be given as second paramater" << endl;
   }


   else // if criterion is number of iteration
      return current_ite<nb_it ? true : false
}

