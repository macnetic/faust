#ifndef __FAUST_CONSTRAINT_REAL_HPP__
#define __FAUST_CONSTRAINT_REAL_HPP__

//#include "faust_constraint_real.h"
#include <cstdlib>
#include "faust_exception.h"

template<typename T>
const char * faust_constraint_real<T>::class_name ="faust_constraint_real<T>::";

template<typename T>
faust_constraint_real<T>::faust_constraint_real() : 
   faust_constraint_generic<T>()
{
   set_default_parameter();
}

template<typename T>
faust_constraint_real<T>::faust_constraint_real(
   const faust_constraint_name& constraint_name_, 
   const int nb_rows_, 
   const int nb_cols_) : 
      faust_constraint_generic<T>(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   set_default_parameter();
}

template<typename T>
faust_constraint_real<T>::faust_constraint_real(
   const faust_constraint_name& constraint_name_,  
   const T default_parameter_,
   const int nb_rows_, 
   const int nb_cols_) : 
      faust_constraint_generic<T>(
         constraint_name_,
         nb_rows_,
         nb_cols_), 
         parameter(default_parameter_)
{
   check_constraint_name();
}

template<typename T>
faust_constraint_real<T>::faust_constraint_real(
   const faust_constraint_real& constraint_) : 
      faust_constraint_generic<T>(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols), 
         parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename T>
void faust_constraint_real<T>::check_constraint_name()const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_NORMCOL:
         break;
      case CONSTRAINT_NAME_NORMLIN:
         break;
      default:
         handleError(class_name,"check_constraint_name : cannot create faust_constraint_real objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename T>
void faust_constraint_real<T>::set_default_parameter()
{

   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_NORMCOL:
         parameter = 0.0;
         break;
      case CONSTRAINT_NAME_NORMLIN:
         parameter = 0.0;
         break;
      default:
         handleError(class_name,"set_default_parameter : cannot create faust_constraint_real objet from an faust_constraint object with this->constraint_name");
         break;
   }
}

template<typename T>
void faust_constraint_real<T>::project(faust_matrix & mat)const
{
	switch (this->constraint_name)
   	{
      		case CONSTRAINT_NAME_NORMCOL:
         		prox_normcol(mat,parameter);
         	break;
      		case CONSTRAINT_NAME_NORMLIN:
         		prox_normlin(mat,parameter);
         	break;
      		default:
         		handleError(class_name,"project : invalid constraint name");
         	break;
   }

}

#endif
