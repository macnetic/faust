#ifndef __FAUST_CONSTRAINT_MAT_HPP__
#define __FAUST_CONSTRAINT_MAT_HPP__

//#include "faust_ConstraintMat.h"
#include <iostream>
#include <cstdlib>
#include "faust_MatDense.h"
#include "faust_exception.h"

template<typename FPP,Device DEVICE>
const char * Faust::ConstraintMat<FPP,DEVICE>::class_name="Faust::ConstraintMat<FPP,DEVICE>::";

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat() :
   Faust::ConstraintGeneric<FPP,DEVICE>()
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraint_name_,
   const faust_unsigned_int nb_rows_,
   const faust_unsigned_int nb_cols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraint_name_,
   const Faust::MatDense<FPP,DEVICE> default_parameter_,
   const faust_unsigned_int nb_rows_,
   const faust_unsigned_int nb_cols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_name_,
         nb_rows_,
         nb_cols_),
         parameter(default_parameter_)
{
   check_constraint_name();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const Faust::ConstraintMat<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols),
         parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         break;
      case CONSTRAINT_NAME_SUPP:
         break;
      default:
         handleError(class_name," cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::set_default_parameter()
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         parameter.setZeros();
         break;
      case CONSTRAINT_NAME_SUPP:
         parameter.setZeros();
         break;
      default:
         handleError(class_name,"set_default_parameter : cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}


template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_CONST:
         mat=parameter;
         break;
      case CONSTRAINT_NAME_SUPP:
         Faust::prox_supp(mat,parameter);
         break;
      default:
         handleError(class_name,"project : invalid constraint_name");
         break;
   }
}



#endif
