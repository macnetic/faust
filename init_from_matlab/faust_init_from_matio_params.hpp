#ifndef __FAUST_INIT_FROM_MATIO_PARAMS_HPP__
#define __FAUST_INIT_FROM_MATIO_PARAMS_HPP__

#include "faust_init_from_matio.h"
#include "faust_init_from_matio_mat.h"
//#include "faust_init_from_matio_params.h"
#include "faust_mat.h"

#include <iostream>
#include <vector>

#include "stopping_criterion.h"
#include "faust_constraint_int.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"

using namespace std;
template<typename T> class faust_constraint_mat;
template<typename T> class faust_constraint_real;

template<typename T>
void init_params_palm_from_matiofile(faust_params_palm<T>& params,const char* fileName, const char* variableName)
{
	
#ifdef __COMPILE_GPU__
   typedef faust_cu_mat<T> faust_matrix ;
#else
   typedef faust_mat<T> faust_matrix ;
#endif
	matvar_t* params_var = faust_matio_read_variable(fileName,variableName);
   
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
	faust_mat<T> data_mat,current_fact;
	vector<faust_matrix > init_facts;
	T init_lambda;	
	
	
	
	
	
	
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
			stopping_criterion<T> stop_cri(niter);
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
			params.init_lambda = (T)init_lambda;
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
				init_facts.push_back(faust_matrix(current_fact));	
			}
			params.init_fact=init_facts;	
		}
		if (strcmp(current_fieldName,"cons")==0)
		{
			vector<const faust_constraint_generic*> consS;

			for (int j=0;j<(current_var->dims[1]);j++)
			{	
				current_cons_var = Mat_VarGetCell(current_var,j);
				add_constraint<T>(consS,current_cons_var);	
			}
			params.cons=consS;
				
		}
	}
}


template<typename T>
void init_params_from_matiofile(faust_params<T>& params, const char* fileName, const char* variableName)
{
	matvar_t* params_var = faust_matio_read_variable(fileName,variableName);
   
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
	faust_mat<T> data_mat,current_fact;
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
			stopping_criterion<T> stop_cri(niter);
			params.stop_crit_2facts=stop_cri;
			//cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"niter2")==0)
		{
			//cout<<"niter2="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion<T> stop_cri(niter);
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
			params.init_lambda = (T) init_lambda;
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
					add_constraint<T>(consS,current_cons_var);				
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



template<typename T>
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
	
	int const_type = getTypeConstraint(name_cons.c_str());
	faust_constraint_name cons_name=getEquivalentConstraint(name_cons.c_str());
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
			consS.push_back(new faust_constraint_int(cons_name,int_parameter,cons_dim1,cons_dim2));
			break;
		
		
		
		// CASE REAL
		case 1 :
		T real_parameter;		
			cons_field_var=Mat_VarGetCell(cons_var,1);
			real_parameter =(T) (((double*) cons_field_var->data))[0];
			consS.push_back(new faust_constraint_real<T>(cons_name,real_parameter,cons_dim1,cons_dim2));
			break;
			
		case 2 :
			faust_mat<T> mat_parameter;
			cons_field_var=Mat_VarGetCell(cons_var,1);
			if ( (cons_dim1 != cons_field_var->dims[0]) || (cons_dim2 != cons_field_var->dims[1]) )
			{
				cerr<<"Error faust_init_from_matio::add_constraint : "<<"mat_parameter of the constraint is invalid"<<endl;
				exit(EXIT_FAILURE); 	
			}
			
			init_mat_from_matvar(mat_parameter,cons_field_var);
			consS.push_back(new faust_constraint_mat<T>(cons_name,mat_parameter,cons_dim1,cons_dim2));
			break;	
	}
	
}


	
template<typename T>	
void Display_params(faust_params<T> & params)
{

	cout<<"NFACTS : "<<params.nb_fact<<endl;
	int nbr_iter_2_fact = 0;
	while(params.stop_crit_2facts.do_continue(nbr_iter_2_fact))
	{
		nbr_iter_2_fact++;
	}
	int nbr_iter_global = 0;
	while(params.stop_crit_global.do_continue(nbr_iter_global))
	{
		nbr_iter_global++;
	}
	cout<<"NBR_ITER_2_FACT : "<<nbr_iter_2_fact << endl;
	cout<<"NBR_ITER_GLOBAL : "<<nbr_iter_global << endl;
	cout<<"VERBOSE : "<<params.isVerbose<<endl;
	cout<<"UPDATEWAY : "<<params.isUpdateWayR2L<<endl;
	cout<<"INIT_LAMBDA : "<<params.init_lambda<<endl;
	cout<<"ISFACTSIDELEFT : "<<params.isFactSideLeft<<endl;
	cout<<"DATA : "<<endl;
	params.data.Display();
	cout<<"INIT_FACTS :"<<endl;
	for (int L=0;L<params.init_fact.size();L++)params.init_fact[L].Display();

	cout<<"CONS : nbr "<< params.cons[0].size()<<endl;

	for (int L=0;L<params.cons[0].size();L++)
	{
		for (int jl=0;jl<params.cons.size();jl++)
		{	//string type_cons;
			//type_cons.resize(0);
			//type_cons=getConstraintType((*params.cons[jl][L]).getConstraintType());
			cout<<"type_cont : "<<params.cons[jl][L]->getType()<<" ";
			cout<<(*params.cons[jl][L]).get_constraint_name();
			cout<<" nb_row :"<<(*params.cons[jl][L]).getRows();
			cout<<" nb_col :"<<(*params.cons[jl][L]).getCols();
			


			if (params.cons[jl][L]->isConstraintParameterInt())
			{	
				faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[jl][L]);
				cout<<" parameter :"<<(*const_int).getParameter()<<endl;
			}
			
			else if (params.cons[jl][L]->isConstraintParameterReal())
			{	
				faust_constraint_real<T>* const_real = (faust_constraint_real<T>*)(params.cons[jl][L]);
				cout<<" parameter :"<<(*const_real).getParameter()<<endl;
			}
			
			else if (params.cons[jl][L]->isConstraintParameterMat())
			{	
				faust_constraint_mat<T>* const_mat = (faust_constraint_mat<T>*)(params.cons[jl][L]);
				cout<<" parameter :"<<endl;
				(*const_mat).getParameter().Display();
			}
			
		}
		cout<<endl<<endl;
	}
}


#endif
