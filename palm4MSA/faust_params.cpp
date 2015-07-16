#include "faust_params.h"
#include "stopping_criterion.h"
#include <iostream>
#include "faust_constraint_int.h"

void faust_params::check_constraint_validity()
{
   
   bool verifSize  =    data.getNbRow()     == cons[0][0]->getRows()
                   && cons[0][0]->getCols() == cons[1][0]->getRows()
                   &&   data.getNbCol()     == cons[1][0]->getCols();

   for (int i=1 ; i<nb_fact-1 ; i++) 
      if (isFactSideLeft)
         verifSize  =  verifSize 
                    && cons[1][i-1]->getRows() == cons[1][i]->getCols()
                    && cons[0][i]->getCols()   == cons[1][i]->getRows()
                    &&    data.getNbRow()      == cons[0][i]->getRows();
      else
         verifSize  =  verifSize 
                    && cons[0][i-1]->getCols() == cons[0][i]->getRows()
                    && cons[0][i]->getCols()   == cons[1][i]->getRows()
                    &&    data.getNbCol()      == cons[1][i]->getCols();


   if (!verifSize)
   {
      std::cerr << "Error in faust_params::check_constraint_validity : Size incompatibility in the constraints" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   for (int i=0 ; i<cons.size() ; i++)  
      if (cons[i].size() != nb_fact-1) 
      {
         std::cerr << "The number of constraints is in conflict with the number of factors" << std::endl;
         exit(EXIT_FAILURE);
      }
 
}

faust_params::faust_params(
         const faust_mat& data_,
         const int nb_fact_,
         const std::vector<std::vector<const faust_constraint_generic*> >& cons_,
         const std::vector<faust_mat>& init_fact_,
         const stopping_criterion& stop_crit_2facts_ /* = stopping_criterion() */,
         const stopping_criterion& stop_crit_global_ /* = stopping_criterion() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const bool isFactSideLeft_ /* = false */,
         const faust_real init_lambda_ /* = 1.0 */,
		 const bool isLambdaComputed_ /* = true */) :
            data(data_), 
            nb_fact(nb_fact_), 
            cons(cons_),
            init_fact(init_fact_),
            stop_crit_2facts(stop_crit_2facts_),
            stop_crit_global(stop_crit_global_),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            isFactSideLeft(isFactSideLeft_),
            init_lambda(init_lambda_),
			isLambdaComputed(isLambdaComputed_)/*,
            nb_rows(data_.getNbRow()),
            nb_cols(data_.getNbCol())*/
{
 check_constraint_validity(); 
}


/*void faust_params::Display() const
{
	std::cout<<"data :  nbRow "<<data.getNbRow()<<" NbCol : "<< data.getNbCol()<<std::endl;
	std::cout<<"nb_facts : "<<nb_fact<<std::endl;
	std::cout<<"Constraints :"<<std::endl;
	std::cout<<"Constraints size:"<<cons.size()<<" "<<cons[0].size()<<std::endl;
	for (int j=0;j<cons.size();j++)
	{	
		std::vector<const faust_constraint_generic*> current_line(cons[j]);
		for (int i=0;i<current_line.size();i++)
		{	std::cout<<j<<" "<<i<<std::endl;
			
			std::cout<<"type : "<<(current_line[i]->get_constraint_name());
			const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(current_line[i]);
			std::cout<<" parameter : "<<const_int->getParameter();
			std::cout<<" DIM1 : "<<cons[j][i]->getRows()<<" DIM2 : "<<cons[j][i]->getCols()<<std::endl;
		}
	}	
	
}*/


faust_params::faust_params() : data(0,0),nb_fact(0),cons(std::vector<std::vector<const faust_constraint_generic*> >()),isFactSideLeft(false),isVerbose(false),isUpdateWayR2L(false),init_fact(std::vector<faust_mat>()),init_lambda(1.0)/*,nb_rows(0),nb_cols(0) */
{}


