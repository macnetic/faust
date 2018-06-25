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
#ifndef __FAUST_INIT_PARAMS_FROM_XML_HPP__
#define __FAUST_INIT_PARAMS_FROM_XML_HPP__


#include <stdlib.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <iostream>
#include <vector>
#include <string>
#include "faust_Transform.h"
#include "faust_MatDense.h"

#include "xml_utils.h"


template<typename FPP,Device DEVICE>

void init_params_from_xml(const char * filename,Faust::Params<FPP,DEVICE> & params)
{

	xmlXPathContextPtr ctxt=get_context(filename);



	std::vector<xmlChar*> contentS;
	// the data matrix is not specified in the xml file

	// m_nbRow
	contentS = get_content(BAD_CAST "/hierarchical_fact/nrow/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file nrow tag must be specified one time");

	}
	params.m_nbRow=xmlXPathCastStringToNumber(contentS[0]);



	contentS = get_content(BAD_CAST "/hierarchical_fact/ncol/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file ncol tag must be specified one time");

	}
	params.m_nbCol=xmlXPathCastStringToNumber(contentS[0]);

	// nb_fact
	contentS = get_content(BAD_CAST "/hierarchical_fact/nb_fact/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file nb_fact tag must be specified one time");

	}
	params.m_nbFact=xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"m_nbFact : "<<params.m_nbFact<<std::endl;

	//niter1 (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter1/text()",ctxt);
	if ((contentS.size() > 1))
		handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file optionnal niter1 tag must be specified at most one time");

	if (contentS.size() == 1)
	{
		Faust::StoppingCriterion<FPP> stop_crit1((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit_2facts = stop_crit1;
	}else
		params.stop_crit_2facts = Faust::StoppingCriterion<FPP>(params.defaultNiter1);

	// std::cout<<"niter1 : "<<params.stop_crit_2facts.get_crit()<<std::endl;

	//niter2 (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter2/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file optionnal niter2 tag must be specified at most one time");
	if (contentS.size() == 1)
	{
		Faust::StoppingCriterion<FPP> stop_crit_global((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit_global = stop_crit_global;
	}else
		params.stop_crit_global = Faust::StoppingCriterion<FPP>(params.defaultNiter2);

	// std::cout<<"niter2 : "<<params.stop_crit_global.get_crit()<<std::endl;



	//updateway (optionnal)
	// std::cout<<"avant get_content updateway"<<std::endl;
	contentS = get_content(BAD_CAST "/hierarchical_fact/updatewayR2L/text()",ctxt);
		// std::cout<<"avant updateway"<<std::endl;
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file :  optionnal updateway tag must be specified at most one time");
	if (contentS.size() == 1)
		params.isUpdateWayR2L = (bool) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.isUpdateWayR2L = params.defaultUpdateWayR2L;	// std::cout<<"updateway : "<<params.isUpdateWayR2L<<std::endl;

	//fact_side
	contentS = get_content(BAD_CAST "/hierarchical_fact/factsideleft/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal factsideleft tag must be specified at most one time");
	if (contentS.size() == 1)
		params.isFactSideLeft = (bool) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.isFactSideLeft = params.defaultFactSideLeft;// std::cout<<"factsideleft : "<<params.isUpdateWayR2L<<std::endl;

	//verbose (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/verbose/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal verbose tag must be specified at most one time");
	if (contentS.size() == 1)
	{
		params.isVerbose = (bool) xmlXPathCastStringToNumber(contentS[0]);
	}
		params.isVerbose = params.defaultVerbosity;
	// std::cout<<"verbose : "<<params.isVerbose<<std::endl;

	//init_lambda
	contentS = get_content(BAD_CAST "/hierarchical_fact/init_lambda/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal init_lambda tag must be specified at most one time");
	if (contentS.size() == 1)
		params.init_lambda = (FPP) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.init_lambda = params.defaultLambda;
	// std::cout<<"init_lambda : "<<params.init_lambda<<std::endl;

	//isConstantStepSize (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/constantstepsize/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal constantstepsize tag must be specified at most one time");
	if (contentS.size() == 1)
		params.isConstantStepSize = (bool) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.isConstantStepSize = params.defaultConstantStepSize;


	//stepsize
	contentS = get_content(BAD_CAST "/hierarchical_fact/stepsize/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal step_size tag must be specified at most one time");
	if (contentS.size() == 1)
		params.step_size = (FPP) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.step_size = params.defaultStepSize;


	//constraints
	std::vector<xmlChar*> contentS_dim1;
	std::vector<xmlChar*> contentS_dim2;
	std::vector<xmlChar*> contentS_type;

	std::vector<string> requestS_begin;

	requestS_begin.push_back("/hierarchical_fact/constraints/row1/constraint/");
	requestS_begin.push_back("/hierarchical_fact/constraints/row2/constraint/");
	vector<vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> > consSS;
	for (int j=0;j<2;j++)
	{

		contentS = get_content(BAD_CAST (requestS_begin[j]+"parameter/text()").c_str(),ctxt);
		contentS_dim1 = get_content(BAD_CAST (requestS_begin[j]+"dim1/text()").c_str(),ctxt);
		contentS_dim2 = get_content(BAD_CAST (requestS_begin[j]+"dim2/text()").c_str(),ctxt);
		contentS_type = get_content(BAD_CAST (requestS_begin[j]+"type/text()").c_str(),ctxt);


		if ((( contentS.size() != (params.m_nbFact-1)) || ((contentS_dim1.size()) != (params.m_nbFact-1))) ||  (( contentS_dim2.size() != (params.m_nbFact-1)) || ((contentS_type.size()) != (params.m_nbFact-1))))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file constraint each row tag must have m_nbFact-1 constraint with each specifying parameter , dim1, dim2,type tags");

		std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> constraintS;
		for (int i=0;i<contentS.size();i++)
		{


			add_constraint<FPP,DEVICE>(constraintS,(char*) contentS_type[i],(char*)contentS[i],(char*)contentS_dim1[i],(char*)contentS_dim2[i]);
		}
		consSS.push_back(constraintS);
	}
	params.cons = consSS;

	// Libération de la mémoire
	xmlXPathFreeContext(ctxt);
	
        params.check_constraint_validity();
}


template<typename FPP,Device DEVICE>
void init_palm_params_from_xml(const char * filename,Faust::ParamsPalm<FPP,DEVICE> & params)
{
	xmlXPathContextPtr ctxt=get_context(filename);



	std::vector<xmlChar*> contentS;
	//data
	// contentS = get_content(BAD_CAST "/palm4MSA/data/text()",ctxt);
	// if ( contentS.size() != 1)
	// {
			// handleError("faust_init_params_from_xml",
			// "init_palm_params_from_xml : invalid file data tag must be specified one time");
			// exit(EXIT_FAILURE);
	// }
	// Faust::MatDense<T> data;
	// data.init_from_file((char *) contentS[0]);
	// params.data=data;

	// m_nbFact
	contentS = get_content(BAD_CAST "/palm4MSA/nb_fact/text()",ctxt);
	if ( contentS.size() != 1)
	{
			handleError("faust_init_params_from_xml",
			"init_palm_params_from_xml : invalid file m_nbFact tag must be specified one time");

	}
	params.nbFact=xmlXPathCastStringToNumber(contentS[0]);

	//niter1 (optionnal)
	contentS = get_content(BAD_CAST "/palm4MSA/niter/text()",ctxt);
	if ((contentS.size() > 1))
		handleError("faust_init_params_from_xml",
			"init_palm_params_from_xml : invalid file optionnal niter1 tag must be specified at most one time");

	if (contentS.size() == 1)
	{
		Faust::StoppingCriterion<FPP> stop_crit1((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit = stop_crit1;
	}else
		params.stop_crit = Faust::StoppingCriterion<FPP>(params.defaultNiter);

	// std::cout<<"niter1 : "<<params.stop_crit_2facts.get_crit()<<std::endl;



	// std::cout<<"niter2 : "<<params.stop_crit_global.get_crit()<<std::endl;



	//updateway (optionnal)
	// std::cout<<"avant get_content updateway"<<std::endl;
	contentS = get_content(BAD_CAST "/palm4MSA/updatewayR2L/text()",ctxt);
		// std::cout<<"avant updateway"<<std::endl;
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_palm_params_from_xml : invalid file :  optionnal updateway tag must be specified at most one time");
	if (contentS.size() == 1)
		params.isUpdateWayR2L = (bool) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.isUpdateWayR2L = params.defaultUpdateWayR2L;	// std::cout<<"updateway : "<<params.isUpdateWayR2L<<std::endl;



	//verbose (optionnal)
	contentS = get_content(BAD_CAST "/palm4MSA/verbose/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal verbose tag must be specified at most one time");
	if (contentS.size() == 1)
	{
		params.isVerbose = (bool) xmlXPathCastStringToNumber(contentS[0]);
	}
		params.isVerbose = params.defaultVerbosity;
	// std::cout<<"verbose : "<<params.isVerbose<<std::endl;

	//init_lambda
	contentS = get_content(BAD_CAST "/palm4MSA/init_lambda/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal init_lambda tag must be specified at most one time");
	if (contentS.size() == 1)
		params.init_lambda = (FPP) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.init_lambda = params.defaultLambda;
		//isConstantStepSize (optionnal)
	contentS = get_content(BAD_CAST "/palm4MSA/constantstepsize/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal constantstepsize tag must be specified at most one time");
	if (contentS.size() == 1)
		params.isConstantStepSize = (bool) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.isConstantStepSize = params.defaultConstantStepSize;


	//stepsize
	contentS = get_content(BAD_CAST "/palm4MSA/stepsize/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal step_size tag must be specified at most one time");
	if (contentS.size() == 1)
		params.step_size = (FPP) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.step_size = params.defaultStepSize;

	//constraints
	std::vector<xmlChar*> contentS_dim1;
	std::vector<xmlChar*> contentS_dim2;
	std::vector<xmlChar*> contentS_type;



		contentS = get_content(BAD_CAST ("/palm4MSA/constraints/constraint/parameter/text()"),ctxt);
		contentS_dim1 = get_content(BAD_CAST ("/palm4MSA/constraints/constraint/dim1/text()"),ctxt);
		contentS_dim2 = get_content(BAD_CAST ("/palm4MSA/constraints/constraint/dim2/text()"),ctxt);
		contentS_type = get_content(BAD_CAST ("/palm4MSA/constraints/constraint/type/text()"),ctxt);


		if ((( contentS.size() != (params.nbFact)) || ((contentS_dim1.size()) != (params.nbFact))) ||  (( contentS_dim2.size() != (params.nbFact)) || ((contentS_type.size()) != (params.nbFact))))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file constraint each row tag must have nbFact-1 constraint with each specifying parameter , dim1, dim2,type tags");

		std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> constraintS;
		for (int i=0;i<contentS.size();i++)
		{


			add_constraint<FPP,DEVICE>(constraintS,(char*) contentS_type[i],(char*)contentS[i],(char*)contentS_dim1[i],(char*)contentS_dim2[i]);
		}
		params.cons = constraintS;

		contentS = get_content(BAD_CAST ("/palm4MSA/init_fact/text()"),ctxt);
		if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal init_fact tag must be specified at most one time");
	if (contentS.size() == 1)
	{
		Faust::Transform<FPP,DEVICE> init_fact_faust;
		init_fact_faust.init_from_file((char*) contentS[0]);
		init_fact_faust.get_facts(params.init_fact);
	}
	else
		params.init_factors();

		// params.check_constraint_validity();
}





template<typename FPP,Device DEVICE>
void add_constraint(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & consS,char* type, char * parameter, char* dim1,char* dim2)
{
	int const_type = get_type_constraint(type);


	int cons_dim1,cons_dim2;
	cons_dim1 =(int) atoi(dim1);
	cons_dim2 =(int) atoi(dim2);
	faust_constraint_name cons_name=get_equivalent_constraint(type);
	switch(const_type)
	{
		// INT CONSTRAINT
		case 0:
		{
			int int_parameter;
			int_parameter =(int) atoi(parameter);
			consS.push_back(new Faust::ConstraintInt<FPP,DEVICE>(cons_name,int_parameter,cons_dim1,cons_dim2));
			break;
		}


		// CASE REAL
		case 1 :
		{
			FPP real_parameter;
			real_parameter =(FPP) atof(parameter);
			consS.push_back(new Faust::ConstraintFPP<FPP,DEVICE>(cons_name,real_parameter,cons_dim1,cons_dim2));
			break;
		}
		case 2 :
		{
			Faust::MatDense<FPP,DEVICE> mat_parameter;
			mat_parameter.init_from_file(parameter);

			if ( (cons_dim1 != mat_parameter.getNbRow()) || (cons_dim2 != mat_parameter.getNbCol()) )
			{
				handleError("faust_init_params_from_xml","::add_constraint : mat_parameter of the constraint is invalid");
				exit(EXIT_FAILURE);
			}
			consS.push_back(new Faust::ConstraintMat<FPP,DEVICE>(cons_name,mat_parameter,cons_dim1,cons_dim2));
			break;
		}
		default :
		{
			std::cerr<<"type : "<<type<<std::endl;
			handleError("faust_init_params_from_xml","init_params_from_xml::add_constraint : invalid constraint type");

		}
	}


}

#endif
