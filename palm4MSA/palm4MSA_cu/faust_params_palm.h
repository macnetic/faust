#ifndef __FAUST_PARAMS_PALM_H__
#define __FAUST_PARAMS_PALM_H__

#include "faust_constant.h"
#include <vector>
#include "stopping_criterion.h"
#include "faust_constraint_generic.h"

#ifdef __COMPILE_GPU__
   #include "faust_mat.h"
#else
   #include "faust_cu_mat.h"
#endif

template<typename T>
class faust_params_palm
{
   private:
#ifdef __COMPILE_GPU__
   typedef faust_cu_mat<T> faust_matrix ;
#else
   typedef faust_mat<T> faust_matrix ;
#endif
      

   public:
      faust_params_palm(
         const faust_matrix& data_,
         const int nb_fact_,
         const std::vector<const faust_constraint_generic*>& cons_,
         const std::vector<faust_matrix >& init_fact_,
         const stopping_criterion<T> & stop_crit_ = stopping_criterion<T>(defaultNiter),
         const bool isVerbose_ = defaultVerbosity ,
         const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
         const T init_lambda_ = defaultLambda,
		 const bool constant_step_size_ = defaultConstantStepSize,
		 const T step_size_ = defaultStepSize);

      void check_constraint_validity();
	   faust_params_palm();	
      ~faust_params_palm(){}


   public:
      // Required members
      faust_matrix data;
      int nb_fact; // number of factors
      std::vector<const faust_constraint_generic*> cons; // vector of constraints

      // Optional members (set to default values if not defined)
      std::vector<faust_matrix > init_fact;
      stopping_criterion<T> stop_crit;
      bool isVerbose;
      bool isUpdateWayR2L;
	  bool isConstantStepSize;
	  T step_size;
      T init_lambda;
		
	  void Display() const;
	  void init_factors();	
	  static const int defaultNiter;
	  static const bool defaultVerbosity;	
	  static const bool defaultUpdateWayR2L;
	  static const T defaultLambda;
	  static const bool defaultConstantStepSize;
	  static const T defaultStepSize;
	  private :
	  static const char *  class_name;

     
      /*const int nb_it;   // number of iterations
      // if isStoppingCriterionError then criterion is error else criterion is number of iteration
      bool  isStoppingCriterionError;
      const T errorThreshold;
      // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than 
      int maxIteration;*/
};

#include "faust_params_palm.hpp"


#endif
