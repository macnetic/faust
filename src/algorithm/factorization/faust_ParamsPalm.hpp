/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
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
const char * Faust::ParamsPalm<FPP,DEVICE>::m_className = "Faust::ParamsPalm<FPP,DEVICE>::";

template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::check_constraint_validity()
{

   bool verifSize  =    data.getNbRow()     == cons[0]->get_rows()
                   &&   data.getNbCol()     == cons[nbFact-1]->get_cols();

	for (int i=0 ; i<nbFact-1; i++)
	{
		verifSize =  verifSize
		&& cons[i]->get_cols() == cons[i+1]->get_rows();
	}
   if (!verifSize)
   {
	   handleError(m_className,"check_constraint_validity : Size incompatibility in the constraints");
   }

	if (init_fact.size() != nbFact)
	{
		handleError(m_className,"check_constraint_validity : conflict between the number of initial factors and nbFact ");
	}
}

template<typename FPP,Device DEVICE>
Faust::ParamsPalm<FPP,DEVICE>::ParamsPalm(
         const Faust::MatDense<FPP,DEVICE>& data_,
         const int nbFact_,
         const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE> *>& cons_,
         const std::vector<Faust::MatDense<FPP,DEVICE> > & init_fact_,
         const Faust::StoppingCriterion<FPP> & stop_crit_ /* = Faust::StoppingCriterion() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const FPP init_lambda_ /* = 1.0 */,
		 const bool constant_step_size_,
		 const FPP step_size_) :
            data(data_),
            nbFact(nbFact_),
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
Faust::ParamsPalm<FPP,DEVICE>::ParamsPalm() : data(0,0),nbFact(0),cons(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>()),isVerbose(defaultVerbosity),isUpdateWayR2L(defaultUpdateWayR2L),init_fact(std::vector<Faust::MatDense<FPP,DEVICE> >()),init_lambda(defaultLambda),isConstantStepSize(defaultConstantStepSize),step_size(defaultStepSize){}

template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::init_factors()
{
	init_fact.resize(nbFact);
   if (!isUpdateWayR2L)
   {
      init_fact[0].resize(cons[0]->get_rows(), cons[0]->get_cols());
      init_fact[0].setZeros();
      for (int i=1 ; i<nbFact ; i++)
      {
         init_fact[i].resize(cons[i]->get_rows(), cons[i]->get_cols());
         init_fact[i].setEyes();
      }
   }
   else
   {
      for (int i=0 ; i<nbFact-1 ; i++)
      {
         init_fact[i].resize(cons[i]->get_rows(), cons[i]->get_cols());
         init_fact[i].setEyes();
      }
      init_fact[nbFact-1].resize(cons[nbFact-1]->get_rows(), cons[nbFact-1]->get_cols());
      init_fact[nbFact-1].setZeros();
   }
}


template<typename FPP,Device DEVICE>
void Faust::ParamsPalm<FPP,DEVICE>::Display() const
{
	std::cout<<"NFACTS : "<<nbFact<<std::endl;
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
			std::cout<<"type_cont : "<<cons[L]->get_type()<<" ";
			std::cout<<(*cons[L]).get_constraint_name();
			std::cout<<" nb_row :"<<(*cons[L]).get_rows();
			std::cout<<" nb_col :"<<(*cons[L]).get_cols();


			if (cons[L]->is_constraint_parameter_int())
			{
				Faust::ConstraintInt<FPP,DEVICE>* const_int = (Faust::ConstraintInt<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<(*const_int).get_parameter()<<std::endl;
			}

			else if (cons[L]->is_constraint_parameter_real())
			{
				Faust::ConstraintFPP<FPP,DEVICE>* const_real = (Faust::ConstraintFPP<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<(*const_real).get_parameter()<<std::endl;
			}

			else if (cons[L]->is_constraint_parameter_mat())
			{
				Faust::ConstraintMat<FPP,DEVICE>* const_mat = (Faust::ConstraintMat<FPP,DEVICE>*)(cons[L]);
				std::cout<<" parameter :"<<std::endl;
				(*const_mat).get_parameter().Display();
			}

		}
		std::cout<<std::endl<<std::endl;
}

#endif
