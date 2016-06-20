#ifndef __FAUST_CONSTRAINT_FPP_HPP__
#define __FAUST_CONSTRAINT_FPP_HPP__


#include <cstdlib>
#include "faust_exception.h"

template<typename FPP,Device DEVICE>
const char * Faust::ConstraintFPP<FPP,DEVICE>::class_name ="Faust::ConstraintFPP<FPP,DEVICE>::";

template<typename FPP,Device DEVICE>
Faust::ConstraintFPP<FPP,DEVICE>::ConstraintFPP() :
   Faust::ConstraintGeneric<FPP,DEVICE>()
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintFPP<FPP,DEVICE>::ConstraintFPP(
   const faust_constraint_name& constraint_name_,
   const int nb_rows_,
   const int nb_cols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintFPP<FPP,DEVICE>::ConstraintFPP(
   const faust_constraint_name& constraint_name_,
   const FPP default_parameter_,
   const int nb_rows_,
   const int nb_cols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_name_,
         nb_rows_,
         nb_cols_),
         parameter(default_parameter_)
{
   check_constraint_name();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintFPP<FPP,DEVICE>::ConstraintFPP(
   const Faust::ConstraintFPP<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols),
         parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename FPP,Device DEVICE>
void Faust::ConstraintFPP<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_NORMCOL:
         break;
      case CONSTRAINT_NAME_NORMLIN:
         break;
      default:
         handleError(class_name,"check_constraint_name : cannot create Faust::ConstraintFPP objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintFPP<FPP,DEVICE>::set_default_parameter()
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
         handleError(class_name,"set_default_parameter : cannot create Faust::ConstraintFPP objet from an faust_constraint object with this->constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintFPP<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat)const
{
	switch (this->constraint_name)
   	{
      		case CONSTRAINT_NAME_NORMCOL:
         		Faust::prox_normcol(mat,parameter);
         	break;
      		case CONSTRAINT_NAME_NORMLIN:
         		Faust::prox_normlin(mat,parameter);
         	break;
      		default:
         		handleError(class_name,"project : invalid constraint name");
         	break;
   }

}

#endif
