#include "faust_ConstraintInt.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"


template<typename FPP,Device DEVICE>
const char * Faust::ConstraintInt<FPP,DEVICE>::class_name="Faust::ConstraintInt";

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt() :
   Faust::ConstraintGeneric<FPP,DEVICE>()
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const faust_constraint_name& constraint_name_,
   const faust_unsigned_int nb_rows_,
   const faust_unsigned_int nb_cols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   std::cout<<nb_rows_<<std::endl;
   std::cout<<nb_cols_<<std::endl;
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const faust_constraint_name& constraint_name_,
   const faust_unsigned_int default_parameter_,
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
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const Faust::ConstraintInt<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols),
         parameter(constraint_.parameter)
{
   check_constraint_name();
}


template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         break;
      case CONSTRAINT_NAME_SPCOL:
         break;
      case CONSTRAINT_NAME_SPLIN:
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         break;
      case CONSTRAINT_NAME_SP_POS:
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         break;
      default:
		handleError(class_name," cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::set_default_parameter()
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         parameter = 0;
         break;
      case CONSTRAINT_NAME_SPCOL:
         parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLIN:
         parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         parameter = 0;
         break;
      case CONSTRAINT_NAME_SP_POS:
         parameter = 0;
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         parameter = 0;
         break;
      default:
		handleError(class_name,"set_default_parameter : cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with this constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         Faust::prox_sp(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPCOL:
         Faust::prox_spcol(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPLIN:
         Faust::prox_splin(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         Faust::prox_splincol(mat,parameter);
         break;
      case CONSTRAINT_NAME_SP_POS:
         Faust::prox_sp_pos(mat,parameter);
         break;
      default:
		handleError(class_name,"project : cannot project with this constraint name");
         break;
   }
}




