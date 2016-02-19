#ifndef __FAUST_CONSTRAINT_MAT_HPP__
#define __FAUST_CONSTRAINT_MAT_HPP__

//#include "faust_constraint_mat.h"
//#include "faust_constraint_mat.h"
#include <iostream>
#include <cstdlib>
#include "faust_mat.h"
#include "faust_exception.h"

template<typename T>
const char * faust_constraint_mat<T>::class_name="faust_constraint_mat<T>::";

template<typename T>
faust_constraint_mat<T>::faust_constraint_mat() : 
   faust_constraint_generic()
{
   set_default_parameter();
}

template<typename T>
faust_constraint_mat<T>::faust_constraint_mat(
   const faust_constraint_name& constraint_name_, 
   const faust_unsigned_int nb_rows_, 
   const faust_unsigned_int nb_cols_) : 
      faust_constraint_generic(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   set_default_parameter();
}

template<typename T>
faust_constraint_mat<T>::faust_constraint_mat(
   const faust_constraint_name& constraint_name_,  
   const faust_matrix default_parameter_,
   const faust_unsigned_int nb_rows_, 
   const faust_unsigned_int nb_cols_) : 
      faust_constraint_generic(
         constraint_name_,
         nb_rows_,
         nb_cols_), 
         parameter(default_parameter_)
{
   check_constraint_name();
}

template<typename T>
faust_constraint_mat<T>::faust_constraint_mat(
   const faust_constraint_mat& constraint_) : 
      faust_constraint_generic(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols), 
         parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename T>
void faust_constraint_mat<T>::check_constraint_name()const
{
   switch (constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         break;
      case CONSTRAINT_NAME_SUPP:
         break;
      default:
         handleError(class_name," cannot create faust_constraint_mat objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename T>
void faust_constraint_mat<T>::set_default_parameter()
{
   switch (constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         parameter.setZeros();
         break;
      case CONSTRAINT_NAME_SUPP:
         parameter.setZeros();
         break;
      default:
         handleError(class_name,"set_default_parameter : cannot create faust_constraint_mat objet from an faust_constraint object with this constraint_name");
         break;
   }
}

#endif
