#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include <iostream>
#include <stdexcept>
#include "faust_exception.h"

template<typename T>
const char * faust_params_palm<T>::class_name = "faust_params_palm<T>::";

template<typename T>
void faust_params_palm<T>::check_constraint_validity()
{
   
   bool verifSize  =    data.getNbRow()     == cons[0]->getRows()
                   &&   data.getNbCol()     == cons[nb_fact-1]->getCols();

	for (int i=0 ; i<nb_fact-1; i++)
	{
		verifSize =  verifSize
		&& cons[i]->getCols() == cons[i+1]->getRows();
	}
   if (!verifSize)
   {
	   handleError(class_name,"check_constraint_validity : Size incompatibility in the constraints");
   }
	
	if (init_fact.size() != nb_fact)
	{
		handleError(class_name,"check_constraint_validity : conflict between the number of initial factors and nb_fact ");
	}
}

template<typename T>
faust_params_palm<T>::faust_params_palm(
         const faust_mat<T>& data_,
         const int nb_fact_,
         const std::vector<const faust_constraint_generic*>& cons_,
         const std::vector<faust_mat<T> > & init_fact_,
         const stopping_criterion<T> & stop_crit_ /* = stopping_criterion() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const T init_lambda_ /* = 1.0 */,
		 const bool isLambdaComputed_ /* true */) :
            data(data_), 
            nb_fact(nb_fact_), 
            cons(cons_),
            init_fact(init_fact_),
            stop_crit(stop_crit_),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            init_lambda(init_lambda_),
			isLambdaComputed(isLambdaComputed_)

{
 check_constraint_validity(); 
}


