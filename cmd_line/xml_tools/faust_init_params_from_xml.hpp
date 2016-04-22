#ifndef __FAUST_INIT_PARAMS_FROM_XML_HPP__
#define __FAUST_INIT_PARAMS_FROM_XML_HPP__


#include <stdlib.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <iostream>
#include <vector>
#include<string>
#include "faust_core.h"

#include "xml_utils.h"


template<typename T>
void init_params_from_xml(const char * filename,faust_params<T> & params)
{

	xmlXPathContextPtr ctxt=get_context(filename);



	std::vector<xmlChar*> contentS;
	// the data matrix is not specified in the xml file
	
	// nb_fact 
	contentS = get_content(BAD_CAST "/hierarchical_fact/nb_fact/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{	
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file nb_fact tag must be specified one time");

	}
	params.nb_fact=xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"nb_fact : "<<params.nb_fact<<std::endl;
	
	//niter1 (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter1/text()",ctxt);
	if ((contentS.size() > 1))
		handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file optionnal niter1 tag must be specified at most one time");
	
	if (contentS.size() == 1)
	{	
		stopping_criterion<T> stop_crit1((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit_2facts = stop_crit1;
	}else
		params.stop_crit_2facts = stopping_criterion<T>(params.defaultNiter1);
	
	// std::cout<<"niter1 : "<<params.stop_crit_2facts.get_crit()<<std::endl;
	
	//niter2 (optionnal)
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter2/text()",ctxt);
	if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file optionnal niter2 tag must be specified at most one time");
	if (contentS.size() == 1)
	{	
		stopping_criterion<T> stop_crit_global((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit_global = stop_crit_global;	
	}else
		params.stop_crit_global = stopping_criterion<T>(params.defaultNiter2);	
		
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
		params.init_lambda = (T) xmlXPathCastStringToNumber(contentS[0]);
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
		params.step_size = (T) xmlXPathCastStringToNumber(contentS[0]);
	else
		params.step_size = params.defaultStepSize;

	
	//constraints
	std::vector<xmlChar*> contentS_dim1;
	std::vector<xmlChar*> contentS_dim2;
	std::vector<xmlChar*> contentS_type;
	
	std::vector<string> requestS_begin;
	
	requestS_begin.push_back("/hierarchical_fact/constraints/row1/constraint/");
	requestS_begin.push_back("/hierarchical_fact/constraints/row2/constraint/");
	vector<vector<const faust_constraint_generic<T>*> > consSS;
	for (int j=0;j<2;j++)
	{	
	
		contentS = get_content(BAD_CAST (requestS_begin[j]+"parameter/text()").c_str(),ctxt);
		contentS_dim1 = get_content(BAD_CAST (requestS_begin[j]+"dim1/text()").c_str(),ctxt);
		contentS_dim2 = get_content(BAD_CAST (requestS_begin[j]+"dim2/text()").c_str(),ctxt);
		contentS_type = get_content(BAD_CAST (requestS_begin[j]+"type/text()").c_str(),ctxt);

	
		if ((( contentS.size() != (params.nb_fact-1)) || ((contentS_dim1.size()) != (params.nb_fact-1))) ||  (( contentS_dim2.size() != (params.nb_fact-1)) || ((contentS_type.size()) != (params.nb_fact-1))))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file constraint each row tag must have nb_fact-1 constraint with each specifying parameter , dim1, dim2,type tags");
		
		std::vector<const faust_constraint_generic<T>*> constraintS;
		for (int i=0;i<contentS.size();i++)
		{		

			
			add_constraint<T>(constraintS,(char*) contentS_type[i],(char*)contentS[i],(char*)contentS_dim1[i],(char*)contentS_dim2[i]);	
		}
		consSS.push_back(constraintS);
	}
	params.cons = consSS;
	
	// Libération de la mémoire
	xmlXPathFreeContext(ctxt);
	// params.check_constraint_validity();
}


template<typename T>
void init_palm_params_from_xml(const char * filename,faust_params_palm<T> & params)
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
	// faust_mat<T> data;
	// data.init_from_file((char *) contentS[0]);
	// params.data=data;
	
	// nb_fact 
	contentS = get_content(BAD_CAST "/palm4MSA/nb_fact/text()",ctxt);
	if ( contentS.size() != 1)
	{	
			handleError("faust_init_params_from_xml",
			"init_palm_params_from_xml : invalid file nb_fact tag must be specified one time");

	}
	params.nb_fact=xmlXPathCastStringToNumber(contentS[0]);
	
	//niter1 (optionnal)
	contentS = get_content(BAD_CAST "/palm4MSA/niter/text()",ctxt);
	if ((contentS.size() > 1))
		handleError("faust_init_params_from_xml",
			"init_palm_params_from_xml : invalid file optionnal niter1 tag must be specified at most one time");
	
	if (contentS.size() == 1)
	{	
		stopping_criterion<T> stop_crit1((int) xmlXPathCastStringToNumber(contentS[0]));
		params.stop_crit = stop_crit1;
	}else
		params.stop_crit = stopping_criterion<T>(params.defaultNiter);
	
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
		params.init_lambda = (T) xmlXPathCastStringToNumber(contentS[0]);
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
		params.step_size = (T) xmlXPathCastStringToNumber(contentS[0]);
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

	
		if ((( contentS.size() != (params.nb_fact)) || ((contentS_dim1.size()) != (params.nb_fact))) ||  (( contentS_dim2.size() != (params.nb_fact)) || ((contentS_type.size()) != (params.nb_fact))))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file constraint each row tag must have nb_fact-1 constraint with each specifying parameter , dim1, dim2,type tags");
		
		std::vector<const faust_constraint_generic<T>*> constraintS;
		for (int i=0;i<contentS.size();i++)
		{		

			
			add_constraint<T>(constraintS,(char*) contentS_type[i],(char*)contentS[i],(char*)contentS_dim1[i],(char*)contentS_dim2[i]);	
		}
		params.cons = constraintS;
		
		contentS = get_content(BAD_CAST ("/palm4MSA/init_fact/text()"),ctxt);
		if ((contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file, optionnal init_fact tag must be specified at most one time");
	if (contentS.size() == 1)
	{		
		faust_core<T> init_fact_faust;
		init_fact_faust.init_from_file((char*) contentS[0]);
		init_fact_faust.get_facts(params.init_fact);
	}
	else
		params.init_factors();
		
		// params.check_constraint_validity(); 
}





template<typename T>
void add_constraint(std::vector<const faust_constraint_generic<T>*> & consS,char* type, char * parameter, char* dim1,char* dim2)
{
	int const_type = getTypeConstraint(type);

	
	int cons_dim1,cons_dim2;
	cons_dim1 =(int) atoi(dim1);
	cons_dim2 =(int) atoi(dim2);
	faust_constraint_name cons_name=getEquivalentConstraint(type);
	switch(const_type)
	{	
		// INT CONSTRAINT
		case 0:
		{
			int int_parameter;		
			int_parameter =(int) atoi(parameter);		
			consS.push_back(new faust_constraint_int<T>(cons_name,int_parameter,cons_dim1,cons_dim2));
			break;
		}
		
		
		// CASE REAL
		case 1 :
		{
			T real_parameter;		
			real_parameter =(T) atof(parameter);
			consS.push_back(new faust_constraint_real<T>(cons_name,real_parameter,cons_dim1,cons_dim2));
			break;
		}	
		case 2 :
		{
			faust_mat<T> mat_parameter;
			mat_parameter.init_from_file(parameter);
			
			if ( (cons_dim1 != mat_parameter.getNbRow()) || (cons_dim2 != mat_parameter.getNbCol()) )
			{
				handleError("faust_init_params_from_xml","::add_constraint : mat_parameter of the constraint is invalid");
				exit(EXIT_FAILURE); 	
			}
			consS.push_back(new faust_constraint_mat<T>(cons_name,mat_parameter,cons_dim1,cons_dim2));
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
