#include "faust_ConstraintInt.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"


template<typename FPP,Device DEVICE>
const char * Faust::ConstraintInt<FPP,DEVICE>::m_className="Faust::ConstraintInt";

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt() :
   Faust::ConstraintGeneric<FPP,DEVICE>()
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const faust_constraint_name& constraintName_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraintName_,
         nbRows_,
         nbCols_)
{
   std::cout<<nbRows_<<std::endl;
   std::cout<<nbCols_<<std::endl;
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const faust_constraint_name& constraintName_,
   const faust_unsigned_int defaultParameter_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraintName_,
         nbRows_,
         nbCols_),
         m_parameter(defaultParameter_)
{
   check_constraint_name();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const Faust::ConstraintInt<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_.constraintName,
         constraint_.nbRows,
         constraint_.nbCols),
         m_parameter(constraint_.parameter)
{
   check_constraint_name();
}


template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->m_constraintName)
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
		handleError(m_className," cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::set_default_parameter()
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPCOL:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLIN:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SP_POS:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         m_parameter = 0;
         break;
      default:
		handleError(m_className,"set_default_parameter : cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with this constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         Faust::prox_sp(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPCOL:
         Faust::prox_spcol(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPLIN:
         Faust::prox_splin(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         Faust::prox_splincol(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SP_POS:
         Faust::prox_sp_pos(mat,m_parameter);
         break;
      default:
		handleError(m_className,"project : cannot project with this constraint name");
         break;
   }
}




