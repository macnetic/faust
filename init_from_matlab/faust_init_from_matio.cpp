#include "faust_init_from_matio.h"
#include "faust_mat.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>
#include "stopping_criterion.h"
#include "faust_constraint_int.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"

using namespace std;
void init_mat_from_matvar(faust_mat & M,matvar_t* var);
void add_constraint(vector<const faust_constraint_generic*> & consS,matvar_t* cons_var);


matvar_t* faust_matio_read_variable(const char* fileName, const char* variableName)
{
   mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDONLY);
   if(matfp == NULL)
   {
      cerr << "error in faust_matio_read_variable : unable to open "<< fileName << endl;
      exit(EXIT_FAILURE);
   }


   matvar_t* matvar = Mat_VarRead(matfp, variableName);
   if(matvar == NULL)
   {
      cerr << "error in faust_matio_read_variable : unable to read "<< variableName <<" from "<< fileName << endl;
      exit(EXIT_FAILURE);
   }
   Mat_Close(matfp);

   return matvar;
}

void init_faust_mat_from_matio_mat(faust_mat& M, const char* fileName, const char* variableName)
{
   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);

   if( matvar->class_type != MAT_C_DOUBLE
       || matvar->rank != 2)
   {
      cerr << "error in init_faust_mat_from_matio_mat : "<< variableName << " seems not to be a matrix." << endl;
      exit(EXIT_FAILURE);
   }
   M.resize(matvar->dims[0], matvar->dims[1]);
   for (size_t i = 0 ; i < matvar->dims[0] * matvar->dims[1] ; i++)
      (M.getData())[i] = (faust_real) (((double*)(matvar->data))[i]);
   Mat_VarFree(matvar);
}







double init_faust_mat_from_matio_double(const char* fileName, const char* variableName)
{
   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);


   if( matvar->class_type != MAT_C_DOUBLE
       || matvar->rank != 2
       || matvar->dims[0] != 1
       || matvar->dims[1] != 1)
   {
      cerr << "error in init_faust_mat_from_matio_double : "<< variableName << " seems not to be a scalar." << endl;
      exit(EXIT_FAILURE);
   }
   double val = ((double*)(matvar->data))[0];
   Mat_VarFree(matvar);
   return val;
}

int init_faust_mat_from_matio_int(const char* fileName, const char* variableName)
{
   double val_d = init_faust_mat_from_matio_double(fileName, variableName);
   int val_i = (int) val_d;
   if ((double)val_i != val_d)
     {
	cerr<<"error in init_faust_mat_from_matio_int : "<< variableName << " seems not to be an integer." <<endl;
	exit(EXIT_FAILURE);
     }
   return val_i;
   
}


bool init_faust_mat_from_matio_bool(const char* fileName, const char* variableName)
{

   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);


   if( (matvar->class_type != MAT_C_DOUBLE && matvar->class_type != MAT_C_UINT8)
       || matvar->rank != 2
       || matvar->dims[0] != 1
       || matvar->dims[1] != 1)
   {
      cerr << "error in init_faust_mat_from_matio_bool : "<< variableName << " seems not to be a scalar." << endl;
      exit(EXIT_FAILURE);
   }
   double val;
   if (matvar->data_type == MAT_T_UINT8)
      val = ((uint8_t*)(matvar->data))[0];
   else if (matvar->data_type == MAT_T_DOUBLE)
      val = ((double*)(matvar->data))[0];
   else
   {
      cerr << "error in init_faust_mat_from_matio_bool : "<< variableName << " wrong data type." << endl;
      exit(EXIT_FAILURE);
   }
   Mat_VarFree(matvar);


   if (val==0.0) 
      return false;
   else if (val == 1)
      return true;
   else
   {
      cerr<<"error in init_faust_mat_from_matio_bool : "<< variableName << " seems not to be a boolean." <<endl;
      exit(EXIT_FAILURE);
   }
}



void write_faust_mat_into_matfile(const faust_mat& M, const char* fileName, const char* variableName)
{
   mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDWR);
   matvar_t *matvar;
   int dim1 = M.getNbRow();
   int dim2 = M.getNbCol();
   double *mat = (double*) malloc(sizeof(double)*dim1*dim2);
   

   if(matfp == NULL)
   {
		matfp = Mat_CreateVer(fileName,NULL,MAT_FT_DEFAULT);
		if ( NULL == matfp ) {
			cerr << "error in write_faust_mat_into_matfile : unable to create "<< fileName << endl;
			 exit(EXIT_FAILURE);
		}
	}
   
   
	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {

        Mat_VarDelete(matfp,matvar->name);
        matvar = NULL;
    }
   

  
	size_t dims[2]={dim1,dim2};
	for (int i = 0 ; i < dim1*dim2; i++) mat[i]=(double)(M.getData())[i];
	
	matvar = Mat_VarCreate(variableName,MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,mat,0);
    if ( NULL == matvar ) {
        cerr<<"error in write_faust_mat_into_matfile : "<< variableName << " unable to create matiovar" <<endl;
		 exit(EXIT_FAILURE);
    } else {
        Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
        Mat_VarFree(matvar);
		
    }
	Mat_Close(matfp);


}




void init_faust_mat_vector_from_matiofile( vector<faust_mat> & vec_M, const char* fileName, const char* variableName)
{
	

	matvar_t* facts_var = faust_matio_read_variable(fileName,"facts");
	cout<<"lecture facts"<<endl;
	matvar_t*   current_fact_var;
	faust_mat current_fact;
	vec_M.resize(0);

	
	if (facts_var->class_type != MAT_C_CELL)
	{
		cerr << "error in init_faust_mat_vector_from_matiofile : facts is'nt a cell_array" << endl;
		exit(EXIT_FAILURE);
	}
	
	int nbr_facts=facts_var->dims[1];
	cout<<"nbr de  facts ="<< nbr_facts<<endl;
	for (int j=0;j<nbr_facts;j++)
	{	

		current_fact_var = Mat_VarGetCell(facts_var,j);
		current_fact.resize(current_fact_var->dims[0],current_fact_var->dims[1]);
		
		for (size_t i = 0 ; i < current_fact_var->dims[0] * current_fact_var->dims[1] ; i++)
      	{
			(((current_fact.getData()))[i]) = (faust_real) (((double*)(current_fact_var->data))[i]);
			
		}
	
		vec_M.push_back(current_fact);	
	}
	
   Mat_VarFree(current_fact_var);
   Mat_VarFree(facts_var);

	
	
	
	
}















void init_params_palm_from_matiofile(faust_params_palm& params,const char* fileName, const char* variableName)
{

	
	matvar_t* params_var = faust_matio_read_variable(fileName,"params");
   
	matvar_t*   current_var;
	matvar_t* current_fact_var;
	matvar_t* current_cons_var;


	
	int nbr_params = Mat_VarGetNumberOfFields(params_var);
	char*const* fieldNames;
	fieldNames = Mat_VarGetStructFieldnames(params_var);
	
	cout<<"nbr de parametre ="<<nbr_params<<endl;
	
	cout<<"***FIELDNAMES*** ="<<endl;
	for (int i=0;i<nbr_params;i++)
	{
		cout<<fieldNames[i]<<endl;
	}
	cout <<endl;
	
	
	char* current_fieldName;
	
	int niter,nfacts,verbose,updateway,dim1,dim2,cons_parameter,cons_dim1,cons_dim2;
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
			cout<<"niter="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit=stop_cri;
			cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"nfacts")==0)
		{
		
			cout<<"nfacts"<<endl;
			nfacts=(int)((double*)(current_var->data))[0];
			cout<<nfacts<<endl;
			params.nb_fact = nfacts;
		}
				
		if (strcmp(current_fieldName,"data")==0)
		{
			cout<<"data"<<endl;
			init_mat_from_matvar(data_mat,current_var);
			data_mat.check_dim_validity();
			data_mat.Display();	
			params.data=data_mat;

		}
		
		if (strcmp(current_fieldName,"verbose")==0)
		{
			cout<<"verbose"<<endl;
			verbose=(int)((double*)(current_var->data))[0];
			cout<<verbose<<endl;
			params.isVerbose = (bool) verbose;
		}
		
		if (strcmp(current_fieldName,"updateway")==0)
		{
			cout<<"updateway"<<endl;
			updateway=(int)((double*)(current_var->data))[0];
			cout<<updateway<<endl;
			params.isUpdateWayR2L=(bool) updateway;
		}
		
		if (strcmp(current_fieldName,"init_lambda")==0)
		{
			cout<<"init_lambda"<<endl;
			init_lambda=(double)((double*)(current_var->data))[0];
			cout<<init_lambda<<endl;
			params.init_lambda = (faust_real)init_lambda;
		}
		if (strcmp(current_fieldName,"init_facts")==0)
		{
			cout<<"init_facts"<<endl;
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
	
	cout<<"nbr de parametre ="<<nbr_params<<endl;
	
	cout<<"***FIELDNAMES*** ="<<endl;
	for (int i=0;i<nbr_params;i++)
	{
		cout<<fieldNames[i]<<endl;
	}
	cout <<endl;
	
	
	char* current_fieldName;
	
	int niter,nfacts,dim1,dim2,cons_parameter,cons_dim1,cons_dim2,fact_side,updateway,verbose;
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
			cout<<"niter1="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit_2facts=stop_cri;
			cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"niter2")==0)
		{
			cout<<"niter2="<<endl;
			niter=(int)((double*)(current_var->data))[0];
			stopping_criterion stop_cri(niter);
			params.stop_crit_global=stop_cri;
			cout<<niter<<endl;
		}
		if (strcmp(current_fieldName,"nfacts")==0)
		{
		
			cout<<"nfacts"<<endl;
			nfacts=(int)((double*)(current_var->data))[0];
			cout<<nfacts<<endl;
			params.nb_fact = nfacts;
		}
				
		if (strcmp(current_fieldName,"data")==0)
		{
			cout<<"data"<<endl;
			init_mat_from_matvar(data_mat,current_var);

			data_mat.Display();
			data_mat.check_dim_validity();		
			params.data=data_mat;

		}
		
		if (strcmp(current_fieldName,"verbose")==0)
		{
			cout<<"verbose"<<endl;
			verbose=(int)((double*)(current_var->data))[0];
			cout<<verbose<<endl;
			params.isVerbose = (bool) verbose;
		}
		
		if (strcmp(current_fieldName,"updateway")==0)
		{
			cout<<"updateway"<<endl;
			updateway=(int)((double*)(current_var->data))[0];
			cout<<updateway<<endl;
			if (updateway != (bool) updateway)
			{     cerr<<"updateway isn't well define"<<endl;
				exit(EXIT_FAILURE);	
			}
			params.isUpdateWayR2L=(bool)updateway;
		}
		
		if (strcmp(current_fieldName,"init_lambda")==0)
		{
			cout<<"init_lambda"<<endl;
			init_lambda=(double)((double*)(current_var->data))[0];
			cout<<init_lambda<<endl;
			params.init_lambda = (faust_real) init_lambda;
			cout<<params.init_lambda<<endl;
		}
		
		if(strcmp(current_fieldName,"fact_side")==0)
		{
			cout<<"fact_side"<<endl;
			fact_side=(int)((double*)(current_var->data))[0];
			cout<<fact_side<<endl;
			params.isFactSideLeft=(bool) fact_side;
		}
		
		if (strcmp(current_fieldName,"cons")==0)
		{
			vector<const faust_constraint_generic*> consS;
			vector<vector<const faust_constraint_generic*> > consSS;
			cout<<"size_tab cont :dim1 "<<current_var->dims[0]<<endl;
			cout<<"size_tab cont :dim2 "<<current_var->dims[1]<<endl;
			int shift = 0;
			if ((current_var->dims[0]) != 2)
			{
				cerr<<"Error init_params_from_matiofile : "<<"params.cons must have 2 rows"<<endl;
				exit(EXIT_FAILURE); 
			}
			
			for (int ll=0;ll<(current_var->dims[0]);ll++)
			{	
				for (int j=0;j<(current_var->dims[1]);j++)
				{	cout<<"j"<<j<<endl;
					current_cons_var = Mat_VarGetCell(current_var,2*j+shift);
					add_constraint(consS,current_cons_var);				
				}
				shift+=1;

				for (int L=0;L<consS.size();L++)
				{
					cout<<consS[L]->get_constraint_name()<<endl;
					//cout<<"params :"<<(*params.cons[L]).getParameter()<<endl;
					cout<<"nb_row :"<<(*consS[L]).getRows()<<endl;
					cout<<"nb_col :"<<(*consS[L]).getCols()<<endl;
					//faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[L]);
					//cout<<"parameter :"<<(*const_int).getParameter()<<endl<<endl;					}
				}
				cout<<endl<<endl;
				
				consSS.push_back(consS);
				consS.resize(0);				
			
			}
			params.cons=consSS;
			
			
				
				
		}
		
	}
	
			
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
	cout<<name_cons<<endl;
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
			int_parameter =(int) (((double*) cons_field_var->data))[0];


					
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


	
	
void init_mat_from_matvar(faust_mat & M,matvar_t* var)
{
	M.resize(var->dims[0],var->dims[1]);
					
	for (size_t k = 0 ; k < var->dims[0] * var->dims[1] ; k++)
	{
		(((M.getData()))[k]) = (faust_real) (((double*)(var->data))[k]);
			
	}
}








