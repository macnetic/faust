#include "faust_constraint_mat.h"
#include "faust_constraint_mat.h"
#include <iostream>
#include <cstdlib>
#include "faust_mat.h"

faust_constraint_mat::faust_constraint_mat() : 
   faust_constraint_generic()
{
   set_default_parameter();
}


faust_constraint_mat::faust_constraint_mat(
   const faust_constraint_name& constraint_name_, 
   const int nb_rows_, 
   const int nb_cols_) : 
      faust_constraint_generic(
         constraint_name_,
         nb_rows_,
         nb_cols_)
{
   set_default_parameter();
}


faust_constraint_mat::faust_constraint_mat(
   const faust_constraint_name& constraint_name_, 
   const int nb_rows_, 
   const int nb_cols_,
   const faust_mat default_parameter_) : 
      faust_constraint_generic(
         constraint_name_,
         nb_rows_,
         nb_cols_), 
         parameter(default_parameter_)
{
   check_constraint_name();
}


faust_constraint_mat::faust_constraint_mat(
   const faust_constraint_mat& constraint_) : 
      faust_constraint_generic(
         constraint_.constraint_name,
         constraint_.nb_rows,
         constraint_.nb_cols), 
         parameter(constraint_.parameter)
{
   check_constraint_name();
}




void faust_constraint_mat::check_constraint_name()const
{
   switch (constraint_name)
   {
      case CONSTRAINT_TYPE_CONST:
         break;
      case CONSTRAINT_TYPE_SUPP:
         break;
      default:
         std::cerr << "Error in faust_constraint_mat::check_constraint_name : cannot create faust_constraint_mat objet from an faust_constraint object with constraint_name= "<< constraint_name << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}

void faust_constraint_mat::set_default_parameter()
{
   switch (constraint_name)
   {
      case CONSTRAINT_TYPE_CONST:
         parameter = 0.0;
         break;
      case CONSTRAINT_TYPE_SUPP:
         parameter = 0.0;
         break;
      default:
         std::cerr << "Error in faust_constraint_mat::set_default_parameter : cannot create faust_constraint_mat objet from an faust_constraint object with constraint_name= "<< constraint_name << std::endl;
         exit(EXIT_FAILURE);
         break;
   }
}
