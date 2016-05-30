#ifndef __FAUST_PARAMS_PALM_HPP__
#define __FAUST_PARAMS_PALM_HPP__

//#include "faust_ParamsPalm.h"
#include "faust_StoppingCriterion.h"
#include <iostream>
#include <stdexcept>
#include "faust_exception.h"

//default value
template<typename FPP,Device DEVICE> const int Faust::ParamsPalm<FPP,DEVICE>::defaultNiter = 300;
template<typename FPP,Device DEVICE> const bool Faust::ParamsPalm<FPP,DEVICE>::defaultVerbosity = false;
template<typename FPP,Device DEVICE> const bool Faust::ParamsPalm<FPP,DEVICE>::defaultUpdateWayR2L = false;
template<typename FPP,Device DEVICE> const FPP Faust::ParamsPalm<FPP,DEVICE>::defaultLambda = 1.0;
template<typename FPP,Device DEVICE> const bool Faust::ParamsPalm<FPP,DEVICE>::defaultConstantStepSize = false;
template<typename FPP,Device DEVICE> const FPP Faust::ParamsPalm<FPP,DEVICE>::defaultStepSize = 1e-16;

template<typename FPP,Device DEVICE>
const char * Faust::ParamsPalm<FPP,DEVICE>::class_name = "Faust::ParamsPalm<FPP,DEVICE>::";

template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::check_constraint_validity()
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

template<typename FPP,Device DEVICE>
Faust::ParamsPalm<FPP,DEVICE>::ParamsPalm(
         const Faust::MatDense<FPP,DEVICE>& data_,
         const int nb_fact_,
         const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE> *>& cons_,
         const std::vector<Faust::MatDense<FPP,DEVICE> > & init_fact_,
         const Faust::StoppingCriterion<FPP> & stop_crit_ /* = Faust::StoppingCriterion() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const FPP init_lambda_ /* = 1.0 */,
		 const bool constant_step_size_,
		 const FPP step_size_) :
            data(data_),
            nb_fact(nb_fact_),
            cons(cons_),
            init_fact(init_fact_),
            stop_crit(stop_crit_),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            init_lambda(init_lambda_),
			isConstantStepSize(constant_step_size_),
			step_size(step_size_)
{
 check_constraint_validity();
}

template<typename FPP,Device DEVICE>
Faust::ParamsPalm<FPP,DEVICE>::ParamsPalm() : data(0,0),nb_fact(0),cons(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>()),isVerbose(defaultVerbosity),isUpdateWayR2L(defaultUpdateWayR2L),init_fact(std::vector<Faust::MatDense<FPP,DEVICE> >()),init_lambda(defaultLambda),isConstantStepSize(defaultConstantStepSize),step_size(defaultStepSize){}

template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::init_factors()
{
	init_fact.resize(nb_fact);
   if (!isUpdateWayR2L)
   {
      init_fact[0].resize(cons[0]->getRows(), cons[0]->getCols());
      init_fact[0].setZeros();
      for (int i=1 ; i<nb_fact ; i++)
      {
         init_fact[i].resize(cons[i]->getRows(), cons[i]->getCols());
         init_fact[i].setEyes();
      }
   }
   else
   {
      for (int i=0 ; i<nb_fact-1 ; i++)
      {
         init_fact[i].resize(cons[i]->getRows(), cons[i]->getCols());
         init_fact[i].setEyes();
      }
      init_fact[nb_fact-1].resize(cons[nb_fact-1]->getRows(), cons[nb_fact-1]->getCols());
      init_fact[nb_fact-1].setZeros();
   }
}


template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::Display() const
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
	std::cout<<"ISCONSTANTSTEPSIZE : "<<isConstantStepSize<<std::endl;
	std::cout<<"step_size : "<<step_size<<std::endl;
	std::cout<<"data :  nbRow "<<data.getNbRow()<<" NbCol : "<< data.getNbCol()<<std::endl;
	std::cout<<"stop_crit : "<<stop_crit.get_crit()<<std::endl;
	/*cout<<"INIT_FACTS :"<<endl;
	for (int L=0;L<init_fact.size();L++)init_fact[L].Display();*/

	std::cout<<"CONSTRAINT  : "<< cons.size()<<std::endl;


		for (unsigned int L=0;L<cons.size();L++)
		{

			//std::string type_cons;
			//type_cons.resize(0);
			//type_cons=getConstraintType((*cons[jl][L]).getConstraintType());
			std::cout<<"type_cont : "<<cons[L]->getType()<<" ";
			std::cout<<(*cons[L]).get_constraint_name();
			std::cout<<" nb_row :"<<(*cons[L]).getRows();
			std::cout<<" nb_col :"<<(*cons[L]).getCols();


			if (cons[L]->isConstraintParameterInt())
			{
				Faust::ConstraintInt<FPP,DEVICE>* const_int = (Faust::ConstraintInt<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<(*const_int).getParameter()<<std::endl;
			}

			else if (cons[L]->isConstraintParameterReal())
			{
				Faust::ConstraintFPP<FPP,DEVICE>* const_real = (Faust::ConstraintFPP<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<(*const_real).getParameter()<<std::endl;
			}

			else if (cons[L]->isConstraintParameterMat())
			{
				Faust::ConstraintMat<FPP,DEVICE>* const_mat = (Faust::ConstraintMat<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<std::endl;
				(*const_mat).getParameter().Display();
			}

		}
		std::cout<<std::endl<<std::endl;
}

#endif
