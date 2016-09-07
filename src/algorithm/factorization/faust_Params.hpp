/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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
#ifndef __FAUST_PARAMS_HPP__
#define __FAUST_PARAMS_HPP__

//#include "Faust::Params.h"
#include "faust_StoppingCriterion.h"
#include <iostream>
#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include <cmath>
#include "faust_exception.h"


template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class ConstraintFPP;
template<typename FPP,Device DEVICE> class ConstraintMat;

template<typename FPP,Device DEVICE>
const char * Faust::Params<FPP,DEVICE>::m_className = "Faust::Params<FPP,DEVICE>::";


/// default Values of Faust::Params
template<typename FPP,Device DEVICE> const bool Faust::Params<FPP,DEVICE>::defaultVerbosity = false;
template<typename FPP,Device DEVICE> const int Faust::Params<FPP,DEVICE>::defaultNiter1 = 500;
template<typename FPP,Device DEVICE> const int Faust::Params<FPP,DEVICE>::defaultNiter2 = 500;
template<typename FPP,Device DEVICE> const bool Faust::Params<FPP,DEVICE>::defaultFactSideLeft = false;
template<typename FPP,Device DEVICE> const bool Faust::Params<FPP,DEVICE>::defaultUpdateWayR2L = false;
template<typename FPP,Device DEVICE> const FPP Faust::Params<FPP,DEVICE>::defaultLambda = 1.0;
template<typename FPP,Device DEVICE> const bool Faust::Params<FPP,DEVICE>::defaultConstantStepSize = false;
template<typename FPP,Device DEVICE> const FPP Faust::Params<FPP,DEVICE>::defaultStepSize = 1e-16;

template<typename FPP,Device DEVICE> const FPP Faust::Params<FPP,DEVICE>::defaultDecreaseSpeed = 1.25;
template<typename FPP,Device DEVICE> const FPP Faust::Params<FPP,DEVICE>::defaultResiduumPercent = 1.4;

template<typename FPP,Device DEVICE>
void Faust::Params<FPP,DEVICE>::check_constraint_validity()
{
    if (cons.size() != 2)
		//handleError("Faust::Params<FPP,DEVICE>::check_constraint_validity :\n cons must have 2 rows instead of %d",cons.size());
        handleError(m_className,"check_constraint_validity :\n cons must have 2 rows");

    for (unsigned int i=0 ; i<cons.size() ; i++)
        if (cons[i].size() != m_nbFact-1)
        {
            //handleError("Faust::Params<FPP,DEVICE>::check_constraint_validity :\n The number of constraints equal to %d is in conflict with the number of factors which is %d\n, number of columns of constraints must be equal to m_nbFact - 1",cons[i].size(),m_nbFact);
            handleError(m_className,"check_constraint_validity :\n The number of constraints equal is in conflict with the number of factors,\n number of columns of constraints must be equal to m_nbFact - 1");
        }

    bool verifSize  =    data.getNbRow()     == cons[0][0]->get_rows()
                    && cons[0][0]->get_cols() == cons[1][0]->get_rows()
                    &&   data.getNbCol()     == cons[1][0]->get_cols();

    for (int i=1 ; i<m_nbFact-1 ; i++)
        if (isFactSideLeft)
            verifSize  =  verifSize
            && cons[1][i-1]->get_rows() == cons[1][i]->get_cols()
            && cons[0][i]->get_cols()   == cons[1][i]->get_rows()
            &&    data.getNbRow()      == cons[0][i]->get_rows();
        else
            verifSize  =  verifSize
            && cons[0][i-1]->get_cols() == cons[0][i]->get_rows()
            && cons[0][i]->get_cols()   == cons[1][i]->get_rows()
            &&    data.getNbCol()      == cons[1][i]->get_cols();


    if (!verifSize)
        handleError(m_className,"Faust::Params<FPP,DEVICE>::check_constraint_validity :\n Size incompatibility in the constraints");

}

template<typename FPP,Device DEVICE>
Faust::Params<FPP,DEVICE>::Params(
	const Faust::MatDense<FPP,DEVICE>& data_,
	const unsigned int nbFact_,
	const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & cons_,
	const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
	const Faust::StoppingCriterion<FPP>& stop_crit_2facts_,
    const Faust::StoppingCriterion<FPP>& stop_crit_global_,
	const FPP residuum_decrease_speed /* = 1.25 */,
	const FPP residuum_prcent /* = 1.4 */,
	const bool isVerbose_ , /* = false */
    const bool isUpdateWayR2L_  , /* = false */
    const bool isFactSideLeft_ , /* = false */
    const FPP init_lambda_  /* = 1.0 */,
	const bool constant_step_size_,
	const FPP step_size_):
        data(data_),
        m_nbFact(nbFact_),
        init_fact(init_fact_),
        stop_crit_2facts(stop_crit_2facts_),
        stop_crit_global(stop_crit_global_),
        isVerbose(isVerbose_),
        isUpdateWayR2L(isUpdateWayR2L_),
        isFactSideLeft(isFactSideLeft_),
        init_lambda(init_lambda_),
        isConstantStepSize(constant_step_size_),
		step_size(step_size_)
{
    if (nbFact_ <= 2)
    {
        //handleError("Faust::Params<FPP,DEVICE>::constructor : 	the number of factor is smaller than 2, use another constructor\n");
        handleError(m_className,"check_constraint_validity : Size incompatibility in the constraints");
    }
    if  (residuum_decrease_speed<=1)
    {
        handleError(m_className,"constructor : residuum_decrease_speed must be strictly greater than  1");
    }
    if ((residuum_prcent<0))
    {
        handleError(m_className,"constructor : residuum_percent must strictly positive");
    }
    if (m_nbFact != cons_.size())
    {
        handleError(m_className,"constructor : m_nbFact and cons_.size() are in conflict\n");
    }

    std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> residuumS_cons;
	std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> factorS_cons;
	double cons_res_parameter = residuum_prcent;
	if(isFactSideLeft)
	{
		for (int i=1;i<m_nbFact-1;i++)
		{
			if (i==1)
			{
                residuumS_cons.push_back(new Faust::ConstraintInt<FPP,DEVICE>(CONSTRAINT_NAME_SP,data.getNbRow()*cons_[m_nbFact-i]->get_rows(),data.getNbRow(),cons_[m_nbFact-i]->get_rows()));
				factorS_cons.push_back(cons_[m_nbFact-i]);
			}else
			{
				std::cout<<m_nbFact-i<<std::endl;
				residuumS_cons.push_back(new Faust::ConstraintInt<FPP,DEVICE>(CONSTRAINT_NAME_SP,std::floor(cons_res_parameter*data.getNbRow()*cons_[m_nbFact-i]->get_rows()+0.5),data.getNbRow(),cons_[m_nbFact-i]->get_rows()));
				std::cout<<m_nbFact-i<<std::endl;
				factorS_cons.push_back(cons_[m_nbFact-i]);
				std::cout<<m_nbFact-i<<std::endl;
			}

			cons_res_parameter=cons_res_parameter/residuum_decrease_speed;
		}
		residuumS_cons.push_back(cons_[0]);

		factorS_cons.push_back(cons_[1]);

		cons.push_back(residuumS_cons);

		cons.push_back(factorS_cons);



	}else
	{
		for (int i=0;i<m_nbFact-2;i++)
		{
			std::cout<<i<<std::endl;
			if (i==0)
			{
				residuumS_cons.push_back(new Faust::ConstraintInt<FPP,DEVICE>(CONSTRAINT_NAME_SP,cons_[i]->get_cols()*data.getNbCol(),cons_[i]->get_cols(),data.getNbCol()));
				factorS_cons.push_back(cons_[0]);
			}else
			{
				residuumS_cons.push_back(new Faust::ConstraintInt<FPP,DEVICE>(CONSTRAINT_NAME_SP,std::floor(cons_res_parameter*cons_[i]->get_cols()*data.getNbCol()+0.5),cons_[i]->get_cols(),data.getNbCol()));
				factorS_cons.push_back(cons_[i]);
			}
				cons_res_parameter=cons_res_parameter/residuum_decrease_speed;

        }

		residuumS_cons.push_back(cons_[m_nbFact-1]);
		factorS_cons.push_back(cons_[m_nbFact-2]);

		cons.push_back(factorS_cons);
		cons.push_back(residuumS_cons);

	}
	check_constraint_validity();

}




template<typename FPP,Device DEVICE>
Faust::Params<FPP,DEVICE>::Params(
         const Faust::MatDense<FPP,DEVICE>& data_,
         const unsigned int nbFact_,
         const std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> >& cons_,
         const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
         const Faust::StoppingCriterion<FPP>& stop_crit_2facts_ /* = Faust::StoppingCriterion<FPP>() */,
         const Faust::StoppingCriterion<FPP>& stop_crit_global_ /* = Faust::StoppingCriterion<FPP>() */,
         const bool isVerbose_ /* = false */,
         const bool isUpdateWayR2L_ /* = false */,
         const bool isFactSideLeft_ /* = false */,
         const FPP init_lambda_ /* = 1.0 */,
		 const bool constant_step_size_ ,
		 const FPP step_size_ ) :
            data(data_),
            m_nbFact(nbFact_),
            cons(cons_),
            init_fact(init_fact_),
            stop_crit_2facts(stop_crit_2facts_),
            stop_crit_global(stop_crit_global_),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            isFactSideLeft(isFactSideLeft_),
            init_lambda(init_lambda_),
			isConstantStepSize(constant_step_size_),
			step_size(step_size_)

{
 check_constraint_validity();
}








template<typename FPP,Device DEVICE>
Faust::Params<FPP,DEVICE>::Params() : data((faust_unsigned_int)0,(faust_unsigned_int)0),m_nbFact(0),cons(std::vector<std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> >()),isFactSideLeft(defaultFactSideLeft),isVerbose(defaultVerbosity),isUpdateWayR2L(defaultUpdateWayR2L),init_fact(std::vector<Faust::MatDense<FPP,DEVICE> >()),init_lambda(defaultLambda),isConstantStepSize(defaultConstantStepSize),step_size(defaultStepSize)
{}


template<typename FPP,Device DEVICE>
void Faust::Params<FPP,DEVICE>::Display() const
{
	std::cout<<"NFACTS : "<<m_nbFact<<std::endl;
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
	std::cout<<"ISCONSTANTSTEPSIZE : "<<isConstantStepSize<<std::endl;
	std::cout<<"step_size : "<<step_size<<std::endl;
	std::cout<<"data :  nbRow "<<data.getNbRow()<<" NbCol : "<< data.getNbCol()<<std::endl;
	std::cout<<"stop_crit_2facts : "<<stop_crit_2facts.get_crit()<<std::endl;
	std::cout<<"stop_crit_global : "<<stop_crit_global.get_crit()<<std::endl;

	/*cout<<"INIT_FACTS :"<<endl;
	for (int L=0;L<init_fact.size();L++)init_fact[L].Display();*/

	std::cout<<"CONSTRAINT  : "<< cons[0].size()<<std::endl;

	for (unsigned int jl=0;jl<cons.size();jl++)
	{

		if (jl == 0)
			if (isFactSideLeft)
				std::cout<<"  RESIDUUMS : "<<std::endl;
			else
				std::cout<<"  FACTORS : "<<std::endl;
		else
			if (isFactSideLeft)
				std::cout<<"  FACTORS : "<<std::endl;
			else
				std::cout<<"  RESIDUUMS : "<<std::endl;

		for (unsigned int L=0;L<cons[0].size();L++)
		{

			//std::string type_cons;
			//type_cons.resize(0);
			//type_cons=get_constraint_type((*cons[jl][L]).get_constraint_type());
			std::cout<<"type_cont : "<<cons[jl][L]->get_type()<<" ";
			std::cout<<(*cons[jl][L]).get_constraint_name();
			std::cout<<" nb_row :"<<(*cons[jl][L]).get_rows();
			std::cout<<" nb_col :"<<(*cons[jl][L]).get_cols();


			if (cons[jl][L]->is_constraint_parameter_int())
			{
				Faust::ConstraintInt<FPP,DEVICE>* const_int = (Faust::ConstraintInt<FPP,DEVICE>*)(cons[jl][L]);
				std::cout<<" parameter :"<<(*const_int).get_parameter()<<std::endl;
			}

			else if (cons[jl][L]->is_constraint_parameter_real())
			{
				Faust::ConstraintFPP<FPP,DEVICE>* const_real = (Faust::ConstraintFPP<FPP,DEVICE>*)(cons[jl][L]);
				std::cout<<" parameter :"<<(*const_real).get_parameter()<<std::endl;
			}

			else if (cons[jl][L]->is_constraint_parameter_mat())
			{
				Faust::ConstraintMat<FPP,DEVICE>* const_mat = (Faust::ConstraintMat<FPP,DEVICE>*)(cons[jl][L]);
				std::cout<<" parameter :"<<std::endl;
				(*const_mat).get_parameter().Display();
			}

		}
		std::cout<<std::endl<<std::endl;
	}

}

template<typename FPP,Device DEVICE>
void Faust::Params<FPP,DEVICE>::init_from_file(const char* filename)
{
	char dataFilename[100];
	int niter1,niter2;

	FILE* fp=fopen(filename,"r");
	if (fp == NULL)
	{
		handleError(m_className,"init_from_file : unable to open file");
	}
	if (feof(fp))
	{
		handleError(m_className,"init_from_file : premature end of file");
	}
	fscanf(fp,"%s\n",dataFilename);
	std::cout<<"dataFilename : "<<dataFilename<<std::endl;
	data.init_from_file(dataFilename);
	std::cout<<"data"<<std::endl;
	// if ((data.getNbCol() > 10) || (data.getNbRow())> 10)
		// data.Display();
	// else
		cout<<"data : nbRow "<<data.getNbRow()<<" nbCol "<<data.getNbCol()<<endl;

	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	fscanf(fp,"%d\n", &m_nbFact);
	std::cout<<"m_nbFact : "<<m_nbFact<<std::endl;
	if (feof(fp))
	{
		handleError(m_className,"init_from_file : premature end of file");
	}
	fscanf(fp,"%d\n", &isVerbose);

	std::cout<<"VERBOSE : "<<isVerbose<<std::endl;
	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	fscanf(fp,"%d\n", &isUpdateWayR2L);
	std::cout<<"UPDATEWAY : "<<isUpdateWayR2L<<std::endl;
	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	if (typeid(FPP)==typeid(double))
	{
		fscanf(fp,"%lf\n", &init_lambda);
	}else
	{
		if (typeid(float)==typeid(float))
		{
			fscanf(fp,"%f\n",&init_lambda);
		}
	}

	std::cout<<"INIT_LAMBDA : "<<init_lambda<<std::endl;
	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	fscanf(fp,"%d\n",&isFactSideLeft);
	std::cout<<"ISFACTSIDELEFT : "<<isFactSideLeft<<std::endl;
	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	fscanf(fp,"%d\n",&niter1);
	std::cout<<"niter1 : "<<niter1<<std::endl;
	Faust::StoppingCriterion<FPP> stopcrit2facts(niter1);
	stop_crit_2facts = stopcrit2facts;
	if (feof(fp))
		handleError(m_className,"init_from_file : premature end of file");
	fscanf(fp,"%d\n",&niter2);
	std::cout<<"niter2 : "<<niter2<<std::endl;
	Faust::StoppingCriterion<FPP> stopcritglobal(niter2);
	stop_crit_global = stopcritglobal;

	vector<const Faust::ConstraintGeneric<FPP,DEVICE> *> consS;
	vector<vector<const Faust::ConstraintGeneric<FPP,DEVICE> *> > consSS;
	for (int i=0;i<2;i++)
	{

		for (int j=0;j<m_nbFact-1;j++)
		{
			if (feof(fp))
				handleError(m_className,"init_from_file : premature end of file");
			char nameCons[100];
			int consDim1,consDim2;
			char cons_parameter[100];
			fscanf(fp,"%s %d %d %s",nameCons,&consDim1,&consDim2,&cons_parameter);
			fscanf(fp,"\n");
			std::cout<<nameCons<<" "<<consDim1<<" "<<consDim2<<" "<<cons_parameter<<std::endl;
			int const_type = get_type_constraint(nameCons);
			faust_constraint_name cons_name=get_equivalent_constraint(nameCons);

			switch(const_type)
			{
				// INT CONSTRAINT
				case 0:
				{
					int int_parameter;
					int_parameter =atoi(cons_parameter);
					consS.push_back(new Faust::ConstraintInt<FPP,DEVICE>(cons_name,int_parameter,consDim1,consDim2));
					break;
				}


				// CASE REAL
				case 1 :
				{
					FPP real_parameter;
					real_parameter=(FPP) atof(cons_parameter);
					consS.push_back(new Faust::ConstraintFPP<FPP,DEVICE>(cons_name,real_parameter,consDim1,consDim2));
					break;
				}
				case 2 :
				{
					Faust::MatDense<FPP,DEVICE> mat_parameter;
					mat_parameter.init_from_file(cons_parameter);

					if ( (consDim1 != mat_parameter.getNbCol()) || (consDim2 != mat_parameter.getNbRow()) )
					{
						handleError(m_className, "init_from_file : invalide dimension of constraint mat_parameter");
					}
					consS.push_back(new Faust::ConstraintMat<FPP,DEVICE>(cons_name,mat_parameter,consDim1,consDim2));
					break;
				}
				default :
				{
					handleError(m_className, "init_from_file : invalid constraint name");
				}
			}
		}
		consSS.push_back(consS);
		consS.resize(0);
		cons = consSS;

	}
	check_bool_validity();
	check_constraint_validity();


}



template<typename FPP,Device DEVICE>
void Faust::Params<FPP,DEVICE>::check_bool_validity()
{
	if (m_nbFact < 1)
		handleError(m_className, "check_bool_validity : m_nbFact must be strcitly greater than 0");
	if ((isVerbose!=0) && (isVerbose!=1))
		handleError(m_className, "check_bool_validity : boolean isVerbose must be equal to 0 or 1");
	if ((isUpdateWayR2L!=0) && (isUpdateWayR2L!=1))
		handleError(m_className, "check_bool_validity : boolean isUpdateWayR2L must be equal to 0 or 1");
	if ((isFactSideLeft!=0) && (isFactSideLeft!=1))
		handleError(m_className, "check_bool_validity : boolean isFactSideLeft must be equal to 0 or 1");
}



#endif
