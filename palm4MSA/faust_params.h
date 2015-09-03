#ifndef __FAUST_PARAMS_H__
#define __FAUST_PARAMS_H__

#include "faust_constant.h"
#include <vector>
#include "faust_mat.h"
#include "stopping_criterion.h"
#include "faust_constraint_generic.h"


class faust_params
{
   public:
	  
	  faust_params(
	  const faust_mat& data_,
	  const int nb_fact_,
	  const std::vector<const faust_constraint_generic*> & cons_,
	  const std::vector<faust_mat>& init_fact_,
	  const stopping_criterion& stop_crit_2facts_ = stopping_criterion(),
	  const stopping_criterion& stop_crit_global_  = stopping_criterion(),
	  const double residuum_decrease_speed = 1.25,
	  const double residuum_prcent = 1.4,
	  const bool isVerbose_ = false ,
	  const bool isUpdateWayR2L_ = false ,
	  const bool isFactSideLeft_ = false ,
	  const faust_real init_lambda_ = 1.0 ,
	  const bool isLambdaComputed_ = true);
	  
	  	

		
      faust_params(
	const faust_mat& data_,
	const int nb_fact_,
	const std::vector<std::vector<const faust_constraint_generic*> >& cons_,
	const std::vector<faust_mat>& init_fact_,
	const stopping_criterion& stop_crit_2facts_ = stopping_criterion(),
	const stopping_criterion& stop_crit_global_  = stopping_criterion(),
	const bool isVerbose_ = false ,
	const bool isUpdateWayR2L_ = false ,
	const bool isFactSideLeft_ = false ,
	const faust_real init_lambda_ = 1.0 ,
	const bool isLambdaComputed_ = true);
		 
	  faust_params();
		
		

      void check_constraint_validity();

      ~faust_params(){}


   public:
      // Required members
      faust_mat data;
      int nb_fact; // number of factors
      std::vector<std::vector<const faust_constraint_generic*> > cons; // vector of constraints
      std::vector<faust_mat> init_fact;

      // Optional members (set to default values if not defined)
      stopping_criterion stop_crit_2facts;
      stopping_criterion stop_crit_global;
      bool isVerbose;
      bool isUpdateWayR2L;
      bool isFactSideLeft;
      faust_real init_lambda;
	  bool isLambdaComputed;
	  void Display() const;

      //const int nb_rows; // number of rows of the first factor
      //const int nb_cols; // number of columns of the last factor
     
      /*const int nb_it;   // number of iterations
      // if isStoppingCriterionError then criterion is error else criterion is number of iteration
      bool  isStoppingCriterionError;
      const faust_real errorThreshold;
      // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than 
      int maxIteration;*/
};

#endif
