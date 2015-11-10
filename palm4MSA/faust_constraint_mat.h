#ifndef __FAUST_CONSTRAINT_MAT_H__
#define __FAUST_CONSTRAINT_MAT_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"
#include "faust_mat.h"


template<typename T> class faust_mat;

//template<typename parameter_type>
template<typename T>
class faust_constraint_mat : public faust_constraint_generic
{
   public:
      faust_constraint_mat(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_, 
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_,  
         const faust_mat<T> parameter_,
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_mat(const faust_constraint_mat& constraint_);

      faust_mat<T> getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;
 
      ~faust_constraint_mat(){};

   private:
      // parameter of constraint
      faust_mat<T> parameter;
	  static const char * class_name;
    
};

#include "faust_constraint_mat.hpp"

#endif
