#include "faust_constraint_generic.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_type.h"
#include <typeinfo>
#include "faust_exception.h"



const char * faust_constraint_generic::class_name="faust_constraint_generic::"; 

const faust_constraint_name faust_constraint_generic::getConstraintType() const 
{
   return constraint_name;
}






const char* faust_constraint_generic::get_constraint_name()const
{
   switch(constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         return "CONSTRAINT_NAME_SP";
      case CONSTRAINT_NAME_SPCOL:
         return "CONSTRAINT_NAME_SPCOL";
      case CONSTRAINT_NAME_SPLIN:
         return "CONSTRAINT_NAME_SPLIN";
      case CONSTRAINT_NAME_NORMCOL:
         return "CONSTRAINT_NAME_NORMCOL";
      case CONSTRAINT_NAME_SPLINCOL:
         return "CONSTRAINT_NAME_SPLINCOL";
      case CONSTRAINT_NAME_L0PEN:
         return "CONSTRAINT_NAME_L0PEN";
      case CONSTRAINT_NAME_L1PEN:
         return "CONSTRAINT_NAME_L1PEN";
      case CONSTRAINT_NAME_CONST:
         return "CONSTRAINT_NAME_CONST";
      case CONSTRAINT_NAME_WAV:
         return "CONSTRAINT_NAME_WAV";
      case CONSTRAINT_NAME_SP_POS:
         return "CONSTRAINT_NAME_SP_POS";
      case CONSTRAINT_NAME_BLKDIAG:
         return "CONSTRAINT_NAME_BLKDIAG";
      case CONSTRAINT_NAME_SPLIN_TEST:
         return "CONSTRAINT_NAME_SPLIN_TEST";
      case CONSTRAINT_NAME_SUPP:
         return "CONSTRAINT_NAME_SUPP";
      case CONSTRAINT_NAME_NORMLIN:
         return "CONSTRAINT_NAME_NORMLIN";
      case CONSTRAINT_NAME_TOEPLITZ:
         return "CONSTRAINT_NAME_TOEPLITZ";
      default:
         return "unknown constraint name";
   }
}

template char const* faust_constraint_generic::getType<double>() const;
template bool faust_constraint_generic::isConstraintParameterInt<double>() const;
template bool faust_constraint_generic::isConstraintParameterReal<double>() const;
template bool faust_constraint_generic::isConstraintParameterMat<double>() const;

template char const* faust_constraint_generic::getType<float>() const;
template bool faust_constraint_generic::isConstraintParameterInt<float>() const;
template bool faust_constraint_generic::isConstraintParameterReal<float>() const;
template bool faust_constraint_generic::isConstraintParameterMat<float>() const;

