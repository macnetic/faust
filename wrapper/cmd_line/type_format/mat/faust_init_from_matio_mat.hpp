#ifndef __FAUST_INIT_FROM_MATIO_MAT_HPP__
#define __FAUST_INIT_FROM_MATIO_MAT_HPP__

#include "faust_init_from_matio.h"
//#include "faust_init_from_matio_mat.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"

#include <iostream>
#include <vector>

using namespace std;

template<typename FPP,Device DEVICE>
void init_faust_mat_from_matio(Faust::MatDense<FPP,DEVICE>& M, const char* fileName, const char* variableName)
{
   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);

   init_mat_from_matvar(M, matvar);

   Mat_VarFree(matvar);
}

template<typename FPP,Device DEVICE>
void init_faust_spmat_from_matio(Faust::MatSparse<FPP,DEVICE>& S, const char* fileName, const char* variableName)
{

   matvar_t* matvar = faust_matio_read_variable(fileName, variableName);

   init_spmat_from_matvar(S, matvar);

   Mat_VarFree(matvar);

}




template<typename FPP,Device DEVICE>
void write_faust_spmat_into_matfile(const Faust::MatSparse<FPP,DEVICE>& M, const char* fileName, const char* variableName)
{
	mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDWR);
   matvar_t *matvar;




   if(matfp == NULL)
   {
		matfp = Mat_CreateVer(fileName,NULL,MAT_FT_DEFAULT);
		if ( NULL == matfp ) {
			cerr << "error in write_faust_mat<FPP,DEVICE>_into_matfile : unable to create "<< fileName << endl;
			 exit(EXIT_FAILURE);
		}
	}


	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
		if (strcmp(matvar->name,variableName) == 0)
		{
			Mat_VarDelete(matfp,matvar->name);
		}
        matvar = NULL;
    }

    write_faust_spmat_list_into_matvar(M,&matvar,variableName);
		Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);

        Mat_VarFree(matvar);


	Mat_Close(matfp);

}

template<typename FPP,Device DEVICE>
void write_faust_spmat_list_into_matvar(const std::vector<Faust::MatSparse<FPP,DEVICE> >& M,matvar_t** matvar, const char* variableName)
{
	vector< Faust::MatDense<FPP,DEVICE> > list_dense_mat;
	list_dense_mat.resize(M.size());
	for (int i=0;i<M.size();i++)
		list_dense_mat[i]=M[i];
	write_faust_mat_list_into_matvar(list_dense_mat,matvar,variableName);

}

template<typename FPP,Device DEVICE>
void write_faust_mat_list_into_matfile(const std::vector< Faust::MatDense<FPP,DEVICE> >& M, const char* fileName, const char* variableName)
{

   mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDWR);
   matvar_t *matvar;




   if(matfp == NULL)
   {
		matfp = Mat_CreateVer(fileName,NULL,MAT_FT_DEFAULT);
		if ( NULL == matfp ) {
			cerr << "error in write_faust_mat<FPP,DEVICE>_into_matfile : unable to create "<< fileName << endl;
			 exit(EXIT_FAILURE);
		}
	}


	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
		if (strcmp(matvar->name,variableName) == 0)
		{
			Mat_VarDelete(matfp,matvar->name);
		}
        matvar = NULL;
    }

    write_faust_mat_list_into_matvar(M,&matvar,variableName);
		Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);

        Mat_VarFree(matvar);


	Mat_Close(matfp);

}


template<typename FPP,Device DEVICE>
void write_faust_mat_list_into_matvar(const std::vector<Faust::MatDense<FPP,DEVICE> >& M,matvar_t** matvar, const char* variableName)
{
	int nbr_mat=M.size();
	size_t dims[2]={(size_t)1,(size_t)nbr_mat};
	(*matvar) = Mat_VarCreate(variableName,MAT_C_CELL,MAT_T_CELL,2,dims,NULL,0);
	if ( NULL == matvar ) {
        cerr<<"error in write_faust_mat_list_into_matfile : "<< variableName << " unable to create matiovar" <<endl;
		 exit(EXIT_FAILURE);
    }
    matvar_t* cell_element;
    for (int i=0;i< nbr_mat;i++)
    {

		write_faust_mat_into_matvar(M[i],&cell_element,"a");
		Mat_VarSetCell((*matvar),i,cell_element);
	}


}




template<typename FPP,Device DEVICE>
void write_faust_mat_into_matvar(const Faust::MatDense<FPP,DEVICE>& M,matvar_t** matvar, const char* variableName)
{

   int dim1 = M.getNbRow();
   int dim2 = M.getNbCol();

	double *mat = (double*) malloc(sizeof(double)*dim1*dim2);
	size_t dims[2]={(size_t)dim1,(size_t)dim2};
	for (int i = 0 ; i < dim1*dim2; i++) mat[i]=(double)(M(i));

	(*matvar) = Mat_VarCreate(variableName,MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,mat,0);
    if ( NULL == matvar ) {
        cerr<<"error in write_faust_mat_into_matvar : "<< variableName << " unable to create matiovar" <<endl;
		 exit(EXIT_FAILURE);
    }
   delete[] mat ;mat=NULL;


}



template<typename FPP,Device DEVICE>
void write_faust_mat_into_matfile(const Faust::MatDense<FPP,DEVICE>& M, const char* fileName, const char* variableName)
{

   mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDWR);
   matvar_t *matvar;




   if(matfp == NULL)
   {
		matfp = Mat_CreateVer(fileName,NULL,MAT_FT_DEFAULT);
		if ( NULL == matfp ) {
			cerr << "error in write_faust_mat<FPP,DEVICE>_into_matfile : unable to create "<< fileName << endl;
			 exit(EXIT_FAILURE);
		}
	}


	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
		if (strcmp(matvar->name,variableName) == 0)
		{
			Mat_VarDelete(matfp,matvar->name);
		}
        matvar = NULL;
    }

		write_faust_mat_into_matvar(M,&matvar,variableName);
		Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);

        Mat_VarFree(matvar);


	Mat_Close(matfp);

}



template<typename FPP,Device DEVICE>
void init_faust_mat_vector_from_matiofile( vector<Faust::MatDense<FPP,DEVICE> > & vec_M, const char* fileName, const char* variableName)
{


	matvar_t* facts_var = faust_matio_read_variable(fileName,"facts");
	//cout<<"lecture facts"<<endl;
	matvar_t*   current_fact_var;
	Faust::MatDense<FPP,DEVICE> current_fact;
	vec_M.resize(0);


	if (facts_var->class_type != MAT_C_CELL)
	{
		cerr << "error in init_faust_mat<FPP,DEVICE>_vector_from_matiofile : facts is'nt a cell_array" << endl;
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
			(((current_fact.getData()))[i]) = (FPP) (((double*)(current_fact_var->data))[i]);

		}

		vec_M.push_back(current_fact);
	}

   Mat_VarFree(current_fact_var);
   Mat_VarFree(facts_var);

}

template<typename FPP,Device DEVICE>
void init_mat_from_matvar(Faust::MatDense<FPP,DEVICE> & M,matvar_t* var)
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
		(((M.getData()))[k]) = (FPP) (((double*)(var->data))[k]);

	}
}

template<typename FPP,Device DEVICE>
void init_spmat_from_matvar(Faust::MatSparse<FPP,DEVICE>& S, matvar_t* var)
{
   mat_sparse_t* mat_sparse = (mat_sparse_t*)var->data;
   if( var->class_type != MAT_C_SPARSE
       || var->rank != 2)
   {
      cerr << "error in init_faust_mat<FPP,DEVICE>_from_matio : "<< "the variable seems not to be a double sparse matrix." << endl;
      exit(EXIT_FAILURE);
   }

   if( var->dims[1] + 1 != mat_sparse->njc
       || mat_sparse->nir < mat_sparse->ndata
       || mat_sparse->jc[var->dims[1]] != mat_sparse->ndata)
   {
      cerr<<"Error in init_faust_spmat<FPP,DEVICE>_from_matio : incorrect dimensions"<<endl;
      exit(EXIT_FAILURE);
   }

   vector<int> rowind(mat_sparse->ndata);
   vector<int> colind(mat_sparse->ndata);
   vector<FPP> values(mat_sparse->ndata);


   for (size_t i = 0 ; i < mat_sparse->ndata ; i++)
      values[i] = (FPP) (((double*)(mat_sparse->data))[i]);

   int cmpt=0;
   for (int i=0 ; i<var->dims[1] ; i++)
      for (int j = mat_sparse->jc[i] ; j < mat_sparse->jc[i + 1] ; j++)
      {
         rowind[cmpt] = mat_sparse->ir[cmpt] ;
         colind[cmpt] = i ;
         cmpt++;
      }
   S = Faust::MatSparse<FPP,DEVICE>(rowind, colind, values, var->dims[0], var->dims[1]);

   if (cmpt != S.getNonZeros())
   {
      cerr<<"Error in init_faust_spmat<FPP,DEVICE>_from_matio : cmpt != nnz : cmpt="<<cmpt<<" ; nnz="<<S.getNonZeros()<<endl;
      exit(EXIT_FAILURE);
   }

}


#endif
