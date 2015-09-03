#include "faust_params.h"
#include "stopping_criterion.h"
#include <iostream>
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"
#include <cmath>

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
   
   for (unsigned int i=0 ; i<cons.size() ; i++)  
      if (cons[i].size() != nb_fact-1) 
      {
         std::cerr << "The number of constraints is in conflict with the number of factors" << std::endl;
         exit(EXIT_FAILURE);
      }
 
}
	faust_params::faust_params(
	  const faust_mat& data_,
	  const int nb_fact_,
	  const std::vector<const faust_constraint_generic*> & cons_,
	  const std::vector<faust_mat>& init_fact_,
	  const stopping_criterion& stop_crit_2facts_,
      const stopping_criterion& stop_crit_global_,
	   const double residuum_decrease_speed /* = 1.25 */,
	  const double residuum_prcent /* = 1.4 */,
	  const bool isVerbose_ , /* = false */
      const bool isUpdateWayR2L_  , /* = false */
      const bool isFactSideLeft_ , /* = false */
      const faust_real init_lambda_ , /* = 1.0 */
      const bool isLambdaComputed_ /* = true */): 
            data(data_), 
            nb_fact(nb_fact_), 
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
  if (nb_fact_ <= 2)
  {
	 std::cerr << "the number of factor is smaller than 2, use another constructor " << std::endl;
         exit(EXIT_FAILURE);  
  }	  
  
  if  (residuum_decrease_speed<=1)
    {
         std::cerr << "residuum_decrease_speed must be strictly greater than  1" << std::endl;
         exit(EXIT_FAILURE);
    }
  	 
   
   if ((residuum_prcent<0))
   {
        std::cerr << "residuum_percent must strictly positive" << std::endl;
        exit(EXIT_FAILURE);
    }
	
	if (nb_fact != cons_.size())
   {
        std::cerr << "nb_fact and cons_.size() are in conflict" << std::endl;
        exit(EXIT_FAILURE);
    }
	
	std::vector<const faust_constraint_generic*> residuumS_cons;
	std::vector<const faust_constraint_generic*> factorS_cons;
	double cons_res_parameter = residuum_prcent; 
	if(isFactSideLeft)
	{	
		for (int i=1;i<nb_fact-1;i++)
		{	
			if (i==1)
			{
					residuumS_cons.push_back(new faust_constraint_int(CONSTRAINT_NAME_SP,data.getNbRow()*cons_[nb_fact-i]->getRows(),data.getNbRow(),cons_[nb_fact-i]->getRows()));
					factorS_cons.push_back(cons_[nb_fact-i]);
					
			}else
			{	
				std::cout<<nb_fact-i<<std::endl;
				residuumS_cons.push_back(new faust_constraint_int(CONSTRAINT_NAME_SP,std::floor(cons_res_parameter*data.getNbRow()*cons_[nb_fact-i]->getRows()+0.5),data.getNbRow(),cons_[nb_fact-i]->getRows()));
				std::cout<<nb_fact-i<<std::endl;
				factorS_cons.push_back(cons_[nb_fact-i]);	
				std::cout<<nb_fact-i<<std::endl;	
			}
			
			cons_res_parameter=cons_res_parameter/residuum_decrease_speed;
		}
		residuumS_cons.push_back(cons_[0]);
		
		factorS_cons.push_back(cons_[1]);
		
		cons.push_back(residuumS_cons);
		
		cons.push_back(factorS_cons);
		 
	//std::cout<<"res size"<< residuumS_cons.size()<<std::endl;
	//std::cout<<"factor size"<< factorS_cons.size()<<std::endl;
	
		
	}else
	{
		for (int i=0;i<nb_fact-2;i++)
		{	
			std::cout<<i<<std::endl;
			if (i==0)
			{
					residuumS_cons.push_back(new faust_constraint_int(CONSTRAINT_NAME_SP,cons_[i]->getCols()*data.getNbCol(),cons_[i]->getCols(),data.getNbCol()));
					factorS_cons.push_back(cons_[0]);
					
			}else
			{	
				residuumS_cons.push_back(new faust_constraint_int(CONSTRAINT_NAME_SP,std::floor(cons_res_parameter*cons_[i]->getCols()*data.getNbCol()+0.5),cons_[i]->getCols(),data.getNbCol()));
				factorS_cons.push_back(cons_[i]);
			}
				cons_res_parameter=cons_res_parameter/residuum_decrease_speed;
				
		}
		
		residuumS_cons.push_back(cons_[nb_fact-1]);
		factorS_cons.push_back(cons_[nb_fact-2]);
		
		cons.push_back(factorS_cons);
		cons.push_back(residuumS_cons);
			
	}
	check_constraint_validity();	

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









faust_params::faust_params() : data(0,0),nb_fact(0),cons(std::vector<std::vector<const faust_constraint_generic*> >()),isFactSideLeft(false),isVerbose(false),isUpdateWayR2L(false),init_fact(std::vector<faust_mat>()),init_lambda(1.0)/*,nb_rows(0),nb_cols(0) */
{}



void faust_params::Display() const
{
	std::cout<<"NFACTS : "<<nb_fact<<std::endl;
	/*int nbr_iter_2_fact = 0;
	while(params.stop_crit_2facts.do_continue(nbr_iter_2_fact))
	{
		nbr_iter_2_fact++;
	}
	int nbr_iter_global = 0;
	while(params.stop_crit_global.do_continue(nbr_iter_global))
	{
		nbr_iter_global++;
	}
	std::cout<<"NBR_ITER_2_FACT : "<<nbr_iter_2_fact << endl;
	std::cout<<"NBR_ITER_GLOBAL : "<<nbr_iter_global << endl;*/
	std::cout<<"VERBOSE : "<<isVerbose<<std::endl;
	std::cout<<"UPDATEWAY : "<<isUpdateWayR2L<<std::endl;
	std::cout<<"INIT_LAMBDA : "<<init_lambda<<std::endl;
	std::cout<<"ISFACTSIDELEFT : "<<isFactSideLeft<<std::endl;
	std::cout<<"data :  nbRow "<<data.getNbRow()<<" NbCol : "<< data.getNbCol()<<std::endl;
	/*cout<<"INIT_FACTS :"<<endl;
	for (int L=0;L<init_fact.size();L++)init_fact[L].Display();*/

	std::cout<<"CONSTRAINT  : "<< cons[0].size()<<std::endl;
	
	for (unsigned int jl=0;jl<cons.size();jl++)
	{
		
		if (jl == 0)
		{
			if (isFactSideLeft)
			{
				std::cout<<"  RESIDUUMS : "<<std::endl; 
			}else
			{
				std::cout<<"  FACTORS : "<<std::endl; 
			}		
		}else
		{
			if (isFactSideLeft)
			{	
				std::cout<<"  FACTORS : "<<std::endl; 
				
			}else
			{
				std::cout<<"  RESIDUUMS : "<<std::endl; 
			}		
			
		}		
		for (unsigned int L=0;L<cons[0].size();L++)
		{
		
			std::string type_cons;
			type_cons.resize(0);
			type_cons=getConstraintType((*cons[jl][L]).getConstraintType());
			std::cout<<"type_cont : "<<type_cons<<" ";
			std::cout<<(*cons[jl][L]).get_constraint_name();
			std::cout<<" nb_row :"<<(*cons[jl][L]).getRows();
			std::cout<<" nb_col :"<<(*cons[jl][L]).getCols();
			
			
			if (strcmp(type_cons.c_str(),"INT") == 0)
			{	
				faust_constraint_int* const_int = (faust_constraint_int*)(cons[jl][L]);
				std::cout<<" parameter :"<<(*const_int).getParameter()<<std::endl;
			}
			
			if (strcmp(type_cons.c_str(),"FAUST_REAL") == 0)
			{	
				faust_constraint_real* const_real = (faust_constraint_real*)(cons[jl][L]);
				std::cout<<" parameter :"<<(*const_real).getParameter()<<std::endl;
			}
			
			if (strcmp(type_cons.c_str(),"FAUST_MAT") == 0)
			{	
				faust_constraint_mat* const_mat = (faust_constraint_mat*)(cons[jl][L]);
				std::cout<<" parameter :"<<std::endl;
				(*const_mat).getParameter().Display();
			}
			
		}
		std::cout<<std::endl<<std::endl;
	}
}

