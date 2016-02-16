#include <iostream>
#include <cstdlib>
#include "faust_cu_mat.h"
#include "faust_exception.h"

template<typename T>
const char * faust_constraint_mat_cu<T>::class_name="faust_constraint_mat_cu<T>::";

template<typename T>
faust_constraint_mat_cu<T>::faust_constraint_mat_cu() : 
   faust_constraint_generic()
{
   set_default_parameter();
}

template<typename T>
faust_constraint_mat_cu<T>::faust_constraint_mat_cu(
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
faust_constraint_mat_cu<T>::faust_constraint_mat_cu(
   const faust_constraint_name& constraint_name_,  
   const faust_cu_mat<T> default_parameter_,
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
faust_constraint_mat_cu<T>::faust_constraint_mat_cu(
   const faust_constraint_mat_cu& constraint_) : 
      faust_constraint_generic(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols), 
         parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename T>
void faust_constraint_mat_cu<T>::check_constraint_name()const
{
   switch (constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         break;
      case CONSTRAINT_NAME_SUPP:
         break;
      default:
         handleError(class_name," cannot create faust_constraint_mat_cu objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename T>
void faust_constraint_mat_cu<T>::set_default_parameter()
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
         handleError(class_name,"set_default_parameter : cannot create faust_constraint_mat_cu objet from an faust_constraint object with this constraint_name");
         break;
   }
}
