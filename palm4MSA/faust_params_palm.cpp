#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include <iostream>

void faust_params_palm::check_constraint_validity()
{
   
   bool verifSize  =    data.getNbRow()     == cons[0]->getRows()
                   && cons[0]->getCols() == cons[1]->getRows()
                   &&   data.getNbCol()     == cons[1]->getCols();


   if (!verifSize)
   {
      std::cerr << "Error in faust_params_palm::check_constraint_validity : Size incompatibility in the constraints" << std::endl;
      exit(EXIT_FAILURE);
   }
 
}

faust_params_palm::faust_params_palm(
         const faust_mat& data_,
         const int nb_fact_,
         const std::vector<const faust_constraint_generic*>& cons_,
         const std::vector<faust_mat>& init_fact_,
         const stopping_criterion& stop_crit_ /* = stopping_criterion() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const faust_real init_lambda_ /* = 1.0 */) :
            data(data_), 
            nb_fact(nb_fact_), 
            cons(cons_),
            init_fact(init_fact_),
            stop_crit(stop_crit_),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            init_lambda(init_lambda_),
            nb_rows(data_.getNbRow()),
            nb_cols(data_.getNbCol())
{
 check_constraint_validity(); 
}


