#include "faust_constraint_int.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"


template<typename FPP>
const char * faust_constraint_int<FPP>::class_name="faust_constraint_int";

template<typename FPP>
faust_constraint_int<FPP>::faust_constraint_int() : 
   faust_constraint_generic<FPP>()
{
   set_default_parameter();
}

template<typename FPP>
faust_constraint_int<FPP>::faust_constraint_int(
   const faust_constraint_name& constraint_name_, 
   const faust_unsigned_int nb_rows_, 
   const faust_unsigned_int nb_cols_) : 
      faust_constraint_generic<FPP>(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   std::cout<<nb_rows_<<std::endl;
   std::cout<<nb_cols_<<std::endl;
   set_default_parameter();
}

template<typename FPP>
faust_constraint_int<FPP>::faust_constraint_int(
   const faust_constraint_name& constraint_name_,  
   const faust_unsigned_int default_parameter_,
   const faust_unsigned_int nb_rows_, 
   const faust_unsigned_int nb_cols_) : 
      faust_constraint_generic<FPP>(
         constraint_name_,
         nb_rows_,
         nb_cols_), 
         parameter(default_parameter_)
{
   check_constraint_name();
}

template<typename FPP>
faust_constraint_int<FPP>::faust_constraint_int(
   const faust_constraint_int& constraint_) : 
      faust_constraint_generic<FPP>(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols), 
         parameter(constraint_.parameter)
{
   check_constraint_name();
}


template<typename FPP>
void faust_constraint_int<FPP>::check_constraint_name()const
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
		handleError(class_name," cannot create faust_constraint_int objet from an faust_constraint object with constraint with constraint_name");
         break;
   }
}

template<typename FPP>
void faust_constraint_int<FPP>::set_default_parameter()
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
		handleError(class_name,"set_default_parameter : cannot create faust_constraint_int objet from an faust_constraint object with constraint with this constraint_name");	
         break;
   }
}

template<typename FPP>
void faust_constraint_int<FPP>::project(faust_matrix & mat) const
{
   switch (this->constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         prox_sp(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPCOL:
         prox_spcol(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPLIN:
         prox_splin(mat,parameter);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         prox_splincol(mat,parameter);
         break;
      case CONSTRAINT_NAME_SP_POS:
         prox_sp_pos(mat,parameter);
         break;
      default:
		handleError(class_name,"project : cannot project with this constraint name");	
         break;
   }
}




