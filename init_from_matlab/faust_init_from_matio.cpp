#include "faust_init_from_matio.h"
#include "faust_mat.h"
#include "faust_constant.h"
#include <iostream>

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

cout<<variableName<<".class_type="<<matvar->class_type<<endl;

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




