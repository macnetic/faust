#include <stdlib.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <iostream>
#include <vector>
#include<string>

#include "xml_utils.h"


template<typename T>
void init_params_from_xml(const char * filename,faust_params<T> & params)
{

	xmlDocPtr doc;
	xmlNodePtr racine;
	 
	// Ouverture du document
	xmlKeepBlanksDefault(0); // Ignore les noeuds texte composant la mise en forme

	 doc = xmlParseFile(filename);
	// Initialisation de l'environnement XPath
	xmlXPathInit();
	// Création du contexte
	xmlXPathContextPtr ctxt = xmlXPathNewContext(doc); // doc est un xmlDocPtr représentant notre catalogue
	if (ctxt == NULL) {
		handleError("faust_init_params_from_xml",
			"init_params_from_xml :Error while the creation of the xpath context ");
	}



	std::vector<xmlChar*> contentS;
	//data
	contentS = get_content(BAD_CAST "/hierarchical_fact/data/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{	
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file data tag must be specified one time");
			exit(EXIT_FAILURE);
	}
	faust_mat<T> data;
	data.init_from_file((char *) contentS[0]);
	params.data=data;
	
	// nb_fact 
	contentS = get_content(BAD_CAST "/hierarchical_fact/nb_fact/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
	{	
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file nb_fact tag must be specified one time");

	}
	params.nb_fact=xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"nb_fact : "<<params.nb_fact<<std::endl;
	
	//niter1
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter1/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file niter1 tag must be specified one time");
	stopping_criterion<T> stop_crit1((int) xmlXPathCastStringToNumber(contentS[0]));		
	params.stop_crit_2facts = stop_crit1;
	// std::cout<<"niter1 : "<<params.stop_crit_2facts.get_crit()<<std::endl;
	
	//niter2
	contentS = get_content(BAD_CAST "/hierarchical_fact/niter2/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file niter2 tag must be specified one time");
	stopping_criterion<T> stop_crit2((int) xmlXPathCastStringToNumber(contentS[0]));		
	params.stop_crit_global = stop_crit2;
	// std::cout<<"niter2 : "<<params.stop_crit_global.get_crit()<<std::endl;
	
	
	
	//updateway
	// std::cout<<"avant get_content updateway"<<std::endl;
	contentS = get_content(BAD_CAST "/hierarchical_fact/updatewayR2L/text()",ctxt);
		// std::cout<<"avant updateway"<<std::endl;
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file updateway tag must be specified one time");
		
	params.isUpdateWayR2L = (bool) xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"updateway : "<<params.isUpdateWayR2L<<std::endl;
		
	//fact_side
	contentS = get_content(BAD_CAST "/hierarchical_fact/factsideleft/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file factsideleft tag must be specified one time");
			
	params.isFactSideLeft = (bool) xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"factsideleft : "<<params.isUpdateWayR2L<<std::endl;

	//verbose
	contentS = get_content(BAD_CAST "/hierarchical_fact/verbose/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid fileverbose tag must be specified one time");
			
	params.isVerbose = (bool) xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"verbose : "<<params.isVerbose<<std::endl;

	//init_lambda
	contentS = get_content(BAD_CAST "/hierarchical_fact/init_lambda/text()",ctxt);
	if (( contentS.size() == 0)||(contentS.size() > 1))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file init_lambda tag must be specified one time");
			
	params.init_lambda = (T) xmlXPathCastStringToNumber(contentS[0]);
	// std::cout<<"init_lambda : "<<params.init_lambda<<std::endl;	
	
	//constraints
	std::vector<xmlChar*> contentS_dim1;
	std::vector<xmlChar*> contentS_dim2;
	std::vector<xmlChar*> contentS_type;
	
	std::vector<string> requestS_begin;
	
	requestS_begin.push_back("/hierarchical_fact/constraints/row1/constraint/");
	requestS_begin.push_back("/hierarchical_fact/constraints/row2/constraint/");
	vector<vector<const faust_constraint_generic*> > consSS;
	for (int j=0;j<2;j++)
	{	
	
		contentS = get_content(BAD_CAST (requestS_begin[j]+"parameter/text()").c_str(),ctxt);
		contentS_dim1 = get_content(BAD_CAST (requestS_begin[j]+"dim1/text()").c_str(),ctxt);
		contentS_dim2 = get_content(BAD_CAST (requestS_begin[j]+"dim2/text()").c_str(),ctxt);
		contentS_type = get_content(BAD_CAST (requestS_begin[j]+"type/text()").c_str(),ctxt);

	
		if ((( contentS.size() != (params.nb_fact-1)) || ((contentS_dim1.size()) != (params.nb_fact-1))) ||  (( contentS_dim2.size() != (params.nb_fact-1)) || ((contentS_type.size()) != (params.nb_fact-1))))
			handleError("faust_init_params_from_xml",
			"init_params_from_xml : invalid file constraint each row tag must have nb_fact-1 constraint with each specifying paramter , dim1, dim2,type tags");
		
		std::vector<const faust_constraint_generic*> constraintS;
		for (int i=0;i<contentS.size();i++)
		{		

			
			add_constraint<T>(constraintS,(char*) contentS_type[i],(char*)contentS[i],(char*)contentS_dim1[i],(char*)contentS_dim2[i]);	
		}
		consSS.push_back(constraintS);
	}
	params.cons = consSS;
	
	// Libération de la mémoire
	params.Display();
	xmlXPathFreeContext(ctxt);
	params.check_constraint_validity();
}

template<typename T>
void add_constraint(std::vector<const faust_constraint_generic*> & consS,char* type, char * parameter, char* dim1,char* dim2)
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
			consS.push_back(new faust_constraint_int(cons_name,int_parameter,cons_dim1,cons_dim2));
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