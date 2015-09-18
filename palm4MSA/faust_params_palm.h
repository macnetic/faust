#ifndef __FAUST_PARAMS_PALM_H__
#define __FAUST_PARAMS_PALM_H__

#include "faust_constant.h"
#include <vector>
#include "faust_mat.h"
#include "stopping_criterion.h"
#include "faust_constraint_generic.h"


class faust_params_palm
{
   public:

      faust_params_palm(
         const faust_mat& data_,
         const int nb_fact_,
         const std::vector<const faust_constraint_generic*>& cons_,
         const std::vector<faust_mat>& init_fact_,
         const stopping_criterion& stop_crit_ = stopping_criterion(),
         const bool isVerbose_ = false ,
         const bool isUpdateWayR2L_ = false ,
         const faust_real init_lambda_ = 1.0 ,
		 const bool isLambdaComputed_ = true);

      void check_constraint_validity();

      ~faust_params_palm(){}


   public:
      // Required members
      faust_mat data;
      int nb_fact; // number of factors
      std::vector<const faust_constraint_generic*> cons; // vector of constraints

      // Optional members (set to default values if not defined)
      std::vector<faust_mat> init_fact;
      stopping_criterion stop_crit;
      bool isVerbose;
      bool isUpdateWayR2L;
      faust_real init_lambda;
	  bool isLambdaComputed;
	  
	  private :
	  static const char *  class_name;


     
      /*const int nb_it;   // number of iterations
      // if isStoppingCriterionError then criterion is error else criterion is number of iteration
      bool  isStoppingCriterionError;
      const faust_real errorThreshold;
      // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than 
      int maxIteration;*/
};

#endif
