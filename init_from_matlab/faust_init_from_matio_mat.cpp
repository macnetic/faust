#include "faust_init_from_matio.h"
#include "faust_init_from_matio_mat.h"
#include "faust_mat.h"
#include "faust_spmat.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>

using namespace std;


void init_faust_mat_from_matio(faust_mat& M, const char* fileName, const char* variableName)
{
   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);

   init_mat_from_matvar(M, matvar);

   Mat_VarFree(matvar);
}

void init_faust_spmat_from_matio(faust_spmat& S, const char* fileName, const char* variableName)
{

   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);
  
   init_spmat_from_matvar(S, matvar);

   Mat_VarFree(matvar);

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
   

  
	size_t dims[2]={(size_t)dim1,(size_t)dim2};
	for (int i = 0 ; i < dim1*dim2; i++) mat[i]=(double)(M(i));
	
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
	//cout<<"lecture facts"<<endl;
	matvar_t*   current_fact_var;
	faust_mat current_fact;
	vec_M.resize(0);

	
	if (facts_var->class_type != MAT_C_CELL)
	{
		cerr << "error in init_faust_mat_vector_from_matiofile : facts is'nt a cell_array" << endl;
		exit(EXIT_FAILURE);
	}
	
	int nbr_facts=facts_var->dims[1];
	//cout<<"nbr de  facts ="<< nbr_facts<<endl;
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
	
	
void init_mat_from_matvar(faust_mat & M,matvar_t* var)
{
	if( var->class_type != MAT_C_DOUBLE
		|| var->rank != 2
		|| var->data_size != sizeof(double) )
	{
		cerr << "error in init_mat_from_matvar : variable seems not to be a double matrix." << endl;
		exit(EXIT_FAILURE);
	}

	M.resize(var->dims[0],var->dims[1]);
					
	for (size_t k = 0 ; k < var->dims[0] * var->dims[1] ; k++)
	{
		(((M.getData()))[k]) = (faust_real) (((double*)(var->data))[k]);
			
	}
}
	
void init_spmat_from_matvar(faust_spmat& S, matvar_t* var)
{
   mat_sparse_t* mat_sparse = (mat_sparse_t*)var->data;
   if( var->class_type != MAT_C_SPARSE
       || var->rank != 2)
   {
      cerr << "error in init_faust_mat_from_matio : "<< "the variable seems not to be a double sparse matrix." << endl;
      exit(EXIT_FAILURE);
   }
   
   if( var->dims[1] + 1 != mat_sparse->njc 
       || mat_sparse->nir < mat_sparse->ndata 
       || mat_sparse->jc[var->dims[1]] != mat_sparse->ndata)
   {
      cerr<<"Error in init_faust_spmat_from_matio : incorrect dimensions"<<endl;
      exit(EXIT_FAILURE);
   }
  
   vector<int> rowind(mat_sparse->ndata);
   vector<int> colind(mat_sparse->ndata);
   vector<faust_real> values(mat_sparse->ndata);


   for (size_t i = 0 ; i < mat_sparse->ndata ; i++)
      values[i] = (faust_real) (((double*)(mat_sparse->data))[i]);

   int cmpt=0;
   for (int i=0 ; i<var->dims[1] ; i++)
      for (int j = mat_sparse->jc[i] ; j < mat_sparse->jc[i + 1] ; j++)
      {
         rowind[cmpt] = mat_sparse->ir[cmpt] ;
         colind[cmpt] = i ;
         cmpt++;
      }
   S = faust_spmat(rowind, colind, values, var->dims[0], var->dims[1]);

   if (cmpt != S.getNonZeros())
   {
      cerr<<"Error in init_faust_spmat_from_matio : cmpt != nnz : cmpt="<<cmpt<<" ; nnz="<<S.getNonZeros()<<endl;
      exit(EXIT_FAILURE);
   }   

}



