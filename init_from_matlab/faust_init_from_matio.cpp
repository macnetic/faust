#include "faust_init_from_matio.h"
#include "faust_mat.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>

using namespace std;

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




