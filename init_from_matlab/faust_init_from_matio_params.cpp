#include "faust_init_from_matio.h"
#include "faust_init_from_matio_mat.h"
#include "faust_init_from_matio_params.h"
#include "faust_mat.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_constraint_int.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"

using namespace std;


void init_params_palm_from_matiofile(faust_params_palm& params,const char* fileName, const char* variableName)
{
	
	matvar_t* params_var = faust_matio_read_variable(fileName,"params");
   
	matvar_t*   current_var;
	matvar_t* current_fact_var;
	matvar_t* current_cons_var;


	
	int nbr_params = Mat_VarGetNumberOfFields(params_var);
	char*const* fieldNames;
	fieldNames = Mat_VarGetStructFieldnames(params_var);
	
	//cout<<"nbr de parametre ="<<nbr_params<<endl;
	
	//cout<<"***FIELDNAMES*** ="<<endl;
	for (int i=0;i<nbr_params;i++)
	{
		//cout<<fieldNames[i]<<endl;
	}
	//cout <<endl;
	
	
	char* current_fieldName;
	
	int niter,nfacts,verbose,update_way,dim1,dim2,cons_parameter,cons_dim1,cons_dim2;
	faust_mat data_mat,current_fact;
	vector<faust_mat> init_facts;
	faust_real init_lambda;	
	
	
	
	
	
	
	for (int i=0;i<nbr_params;i++)
	{
		current_fieldName = fieldNames[i];	
		int id_champ = -1;

		

		
	
	    
		current_var=Mat_VarGetStructFieldByName(params_var,current_fieldName,0);
		if (current_var == NULL)
		{
			cerr<<"error cannot access to the field : "<<current_fieldName<<endl;
			exit(EXIT_FAILURE);
			
		}
		if (strcmp(current_fieldName,"niter")==0)
		{
			//cout<<"niter="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit=stop_cri;
			//cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"nfacts")==0)
		{
		
			//cout<<"nfacts"<<endl;
			nfacts=(int)((double*)(current_var->data))[0];
			//cout<<nfacts<<endl;
			params.nb_fact = nfacts;
		}
				
		if (strcmp(current_fieldName,"data")==0)
		{
			//cout<<"data"<<endl;
			init_mat_from_matvar(data_mat,current_var);
			data_mat.check_dim_validity();
			data_mat.Display();	
			params.data=data_mat;

		}
		
		if (strcmp(current_fieldName,"verbose")==0)
		{
			//cout<<"verbose"<<endl;
			verbose=(int)((double*)(current_var->data))[0];
			//cout<<verbose<<endl;
			params.isVerbose = (bool) verbose;
		}
		
		if (strcmp(current_fieldName,"update_way")==0)
		{
			//cout<<"update_way"<<endl;
			update_way=(int)((double*)(current_var->data))[0];
			//cout<<update_way<<endl;
			params.isUpdateWayR2L=(bool) update_way;
		}
		
		if (strcmp(current_fieldName,"init_lambda")==0)
		{
			//cout<<"init_lambda"<<endl;
			init_lambda=(double)((double*)(current_var->data))[0];
			//cout<<init_lambda<<endl;
			params.init_lambda = (faust_real)init_lambda;
		}
		if (strcmp(current_fieldName,"init_facts")==0)
		{
			//cout<<"init_facts"<<endl;
			init_facts.resize(0);
			
			for (int j=0;j<(current_var->dims[1]);j++)
			{	

				current_fact_var = Mat_VarGetCell(current_var,j);
				
				init_mat_from_matvar(current_fact,current_fact_var);
				current_fact.check_dim_validity();		
				current_fact.Display();
				init_facts.push_back(current_fact);	
			}
			params.init_fact=init_facts;	
		}
		if (strcmp(current_fieldName,"cons")==0)
		{
			vector<const faust_constraint_generic*> consS;

			for (int j=0;j<(current_var->dims[1]);j++)
			{	
				current_cons_var = Mat_VarGetCell(current_var,j);
				add_constraint(consS,current_cons_var);	
			}
			params.cons=consS;
				
		}
	}
}



void init_params_from_matiofile(faust_params& params, const char* fileName, const char* variableName)
{
	matvar_t* params_var = faust_matio_read_variable(fileName,"params");
   
	matvar_t*   current_var;
	matvar_t* current_cons_var;

	
	int nbr_params = Mat_VarGetNumberOfFields(params_var);
	char*const* fieldNames;
	fieldNames = Mat_VarGetStructFieldnames(params_var);
	
	//cout<<"nbr de parametre ="<<nbr_params<<endl;
	
	//cout<<"***FIELDNAMES*** ="<<endl;
	for (int i=0;i<nbr_params;i++)
	{
		//cout<<fieldNames[i]<<endl;
	}
	//cout <<endl;
	
	
	char* current_fieldName;
	
	int niter,nfacts,dim1,dim2,cons_parameter,cons_dim1,cons_dim2,fact_side,update_way,verbose;
	faust_mat data_mat,current_fact;
	double init_lambda;	
	
	
	
	
	
	
	for (int i=0;i<nbr_params;i++)
	{
		current_fieldName = fieldNames[i];	
		int id_champ = -1;

		

		
	
	    
		current_var=Mat_VarGetStructFieldByName(params_var,current_fieldName,0);
		if (current_var == NULL)
		{
			cerr<<"error cannot access to the field : "<<current_fieldName<<endl;
			exit(EXIT_FAILURE);	
			
		}
		if (strcmp(current_fieldName,"niter1")==0)
		{
			//cout<<"niter1="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit_2facts=stop_cri;
			//cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"niter2")==0)
		{
			//cout<<"niter2="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit_global=stop_cri;
			//cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"nfacts")==0)
		{
		
			//cout<<"nfacts"<<endl;
			nfacts=(int)((double*)(current_var->data))[0];
			//cout<<nfacts<<endl;
			params.nb_fact = nfacts;
		}
				
		if (strcmp(current_fieldName,"data")==0)
		{
			//cout<<"data"<<endl;
			init_mat_from_matvar(data_mat,current_var);

			//data_mat.Display();
			data_mat.check_dim_validity();		
			params.data=data_mat;

		}
		
		if (strcmp(current_fieldName,"verbose")==0)
		{
			//cout<<"verbose"<<endl;
			verbose=(int)((double*)(current_var->data))[0];
			//cout<<verbose<<endl;
			params.isVerbose = (bool) verbose;
		}
		
		if (strcmp(current_fieldName,"update_way")==0)
		{
			//cout<<"update_way"<<endl;
			update_way=(int)((double*)(current_var->data))[0];
			//cout<<update_way<<endl;
			if (update_way != (bool) update_way)
			{     cerr<<"update_way isn't well define"<<endl;
				exit(EXIT_FAILURE);	
			}
			params.isUpdateWayR2L=(bool)update_way;
		}
		
		if (strcmp(current_fieldName,"init_lambda")==0)
		{
			//cout<<"init_lambda"<<endl;
			init_lambda=(double)((double*)(current_var->data))[0];
			//cout<<init_lambda<<endl;
			params.init_lambda = (faust_real) init_lambda;
			//cout<<params.init_lambda<<endl;
		}
		
		if(strcmp(current_fieldName,"fact_side")==0)
		{
			////cout<<"fact_side"<<endl;
			fact_side=(int)((double*)(current_var->data))[0];
			////cout<<fact_side<<endl;
			params.isFactSideLeft=(bool) fact_side;
		}
		
		if (strcmp(current_fieldName,"cons")==0)
		{
			vector<const faust_constraint_generic*> consS;
			vector<vector<const faust_constraint_generic*> > consSS;
			////cout<<"size_tab cont :dim1 "<<current_var->dims[0]<<endl;
			////cout<<"size_tab cont :dim2 "<<current_var->dims[1]<<endl;
			int shift = 0;
			if ((current_var->dims[0]) != 2)
			{
				cerr<<"Error init_params_from_matiofile : "<<"params.cons must have 2 rows"<<endl;
				exit(EXIT_FAILURE); 
			}
			
			for (int ll=0;ll<(current_var->dims[0]);ll++)
			{	
				for (int j=0;j<(current_var->dims[1]);j++)
				{	////cout<<"j"<<j<<endl;
					current_cons_var = Mat_VarGetCell(current_var,2*j+shift);
					add_constraint(consS,current_cons_var);				
				}
				shift+=1;

				for (int L=0;L<consS.size();L++)
				{
					////cout<<consS[L]->get_constraint_name()<<endl;
					////cout<<"params :"<<(*params.cons[L]).getParameter()<<endl;
					////cout<<"nb_row :"<<(*consS[L]).getRows()<<endl;
					////cout<<"nb_col :"<<(*consS[L]).getCols()<<endl;
					//faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[L]);
					////cout<<"parameter :"<<(*const_int).getParameter()<<endl<<endl;
				}
				//cout<<endl<<endl;
				
				consSS.push_back(consS);
				consS.resize(0);				
			
			}
			params.cons=consSS;
			
			
				
				
		}
		
	}
	params.check_constraint_validity();
	
			
}




void add_constraint(std::vector<const faust_constraint_generic*> & consS,matvar_t* cons_var)
{
	
	matvar_t* cons_name_var, *cons_field_var;
	string name_cons;
	cons_name_var = Mat_VarGetCell(cons_var,0);
	name_cons.resize(0);
	for(int k=0;k<cons_name_var->dims[1];k++)
	{
		name_cons+= (char) (((char*)(cons_name_var->data))[k]);
	}
	//cout<<name_cons<<endl;
	bool is_const_int =((strcmp(name_cons.c_str(),"sp") == 0) || (strcmp(name_cons.c_str(),"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"lOpen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"l1pen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"wav") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"blkdiag") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"splin_test") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"supp") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(name_cons.c_str(),"normlin") == 0)));
	
	bool is_const_real = ((strcmp(name_cons.c_str(),"normcol") == 0) || (strcmp(name_cons.c_str(),"normlin")==0));
	bool is_const_mat =  ((strcmp(name_cons.c_str(),"supp") == 0) || (strcmp(name_cons.c_str(),"const")==0));
					
	int const_type = -1;
	if (is_const_int)
	{
		const_type = 0;
	}
	if (is_const_real)
	{
		const_type = 1;
	}
	if (is_const_mat)
	{
		const_type = 2;
	}
					int nbr_field = (cons_var->dims[1]);
				

	if (nbr_field != 4)
	{
		cerr<<"Error faust_init_from_matio::add_constraint : "<<"params.cons{i} must have 4 fields"<<endl;
		exit(EXIT_FAILURE); 
	}
	int cons_dim1,cons_dim2;
	cons_field_var=Mat_VarGetCell(cons_var,2);
	cons_dim1 =(int) (((double*) cons_field_var->data))[0];
	cons_field_var=Mat_VarGetCell(cons_var,3);
	cons_dim2 =(int) (((double*) cons_field_var->data))[0];
	
	switch(const_type)
	{	
		// INT CONSTRAINT
		case 0:

			int int_parameter;		
			cons_field_var=Mat_VarGetCell(cons_var,1);
			int_parameter =(int) round((((double*) cons_field_var->data))[0]);


					
			faust_constraint_name cons_name;			
			if (strcmp(name_cons.c_str(),"sp") == 0)
			{
				cons_name = CONSTRAINT_NAME_SP;
			}

			if (strcmp(name_cons.c_str(),"sppos") == 0)
			{
				cons_name = CONSTRAINT_NAME_SP_POS;
			}

			if (strcmp(name_cons.c_str(),"spcol") == 0)
			{
				cons_name = CONSTRAINT_NAME_SPCOL;
			}			
					
			if (strcmp(name_cons.c_str(),"splin") == 0)
			{
				cons_name = CONSTRAINT_NAME_SPLIN;
			}	
					
			consS.push_back(new faust_constraint_int(cons_name,int_parameter,cons_dim1,cons_dim2));
			break;
		
		
		
		// CASE REAL
		case 1 :
		faust_real real_parameter;		
			cons_field_var=Mat_VarGetCell(cons_var,1);
			real_parameter =(faust_real) (((double*) cons_field_var->data))[0];
			
			if (strcmp(name_cons.c_str(),"normcol") == 0)
			{
				cons_name = CONSTRAINT_NAME_NORMCOL;
			}

			if (strcmp(name_cons.c_str(),"normlin") == 0)
			{
				cons_name = CONSTRAINT_NAME_NORMLIN;
			}
			consS.push_back(new faust_constraint_real(cons_name,real_parameter,cons_dim1,cons_dim2));
			break;
			
		case 2 :
			faust_mat mat_parameter;
			cons_field_var=Mat_VarGetCell(cons_var,1);
			if ( (cons_dim1 != cons_field_var->dims[0]) || (cons_dim2 != cons_field_var->dims[1]) )
			{
				cerr<<"Error faust_init_from_matio::add_constraint : "<<"mat_parameter of the constraint is invalid"<<endl;
				exit(EXIT_FAILURE); 	
			}
			
			init_mat_from_matvar(mat_parameter,cons_field_var);
			
			if (strcmp(name_cons.c_str(),"const") == 0)
			{
				cons_name = CONSTRAINT_NAME_CONST;
			}

			if (strcmp(name_cons.c_str(),"supp") == 0)
			{
				cons_name = CONSTRAINT_NAME_SUPP;
			}
			consS.push_back(new faust_constraint_mat(cons_name,mat_parameter,cons_dim1,cons_dim2));
			break;
			
		/*default :
			cerr<<"Error faust_init_from_matio::add_constraint : "<<"constraint's name "<<name_cons<<" is invalid "<<endl;
			exit(EXIT_FAILURE);
		*/		
	}
	
}


	
	
// void Display_params(faust_params & params)
// {

	// cout<<"NFACTS : "<<params.nb_fact<<endl;
	// int nbr_iter_2_fact = 0;
	// while(params.stop_crit_2facts.do_continue(nbr_iter_2_fact))
	// {
		// nbr_iter_2_fact++;
	// }
	// int nbr_iter_global = 0;
	// while(params.stop_crit_global.do_continue(nbr_iter_global))
	// {
		// nbr_iter_global++;
	// }
	// cout<<"NBR_ITER_2_FACT : "<<nbr_iter_2_fact << endl;
	// cout<<"NBR_ITER_GLOBAL : "<<nbr_iter_global << endl;
	// cout<<"VERBOSE : "<<params.isVerbose<<endl;
	// cout<<"UPDATEWAY : "<<params.isUpdateWayR2L<<endl;
	// cout<<"INIT_LAMBDA : "<<params.init_lambda<<endl;
	// cout<<"ISFACTSIDELEFT : "<<params.isFactSideLeft<<endl;
	// cout<<"DATA : "<<endl;
	// params.data.Display();
	// cout<<"INIT_FACTS :"<<endl;
	// for (int L=0;L<params.init_fact.size();L++)params.init_fact[L].Display();

	// cout<<"CONS : nbr "<< params.cons[0].size()<<endl;

	// for (int L=0;L<params.cons[0].size();L++)
	// {
		// for (int jl=0;jl<params.cons.size();jl++)
		// {	//string type_cons;
			type_cons.resize(0);
			type_cons=getConstraintType((*params.cons[jl][L]).getConstraintType());
			// cout<<"type_cont : "<<params.cons[jl][L]->getType()<<" ";
			// cout<<(*params.cons[jl][L]).get_constraint_name();
			// cout<<" nb_row :"<<(*params.cons[jl][L]).getRows();
			// cout<<" nb_col :"<<(*params.cons[jl][L]).getCols();
			


			// if (params.cons[jl][L]->isConstraintParameterInt())
			// {	
				// faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[jl][L]);
				// cout<<" parameter :"<<(*const_int).getParameter()<<endl;
			// }
			
			// else if (params.cons[jl][L]->isConstraintParameterReal())
			// {	
				// faust_constraint_real* const_real = (faust_constraint_real*)(params.cons[jl][L]);
				// cout<<" parameter :"<<(*const_real).getParameter()<<endl;
			// }
			
			// else if (params.cons[jl][L]->isConstraintParameterMat())
			// {	
				// faust_constraint_mat* const_mat = (faust_constraint_mat*)(params.cons[jl][L]);
				// cout<<" parameter :"<<endl;
				// (*const_mat).getParameter().Display();
			// }
			
		// }
		// cout<<endl<<endl;
	// }
// }






