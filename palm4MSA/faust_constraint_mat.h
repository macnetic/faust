#ifndef __FAUST_CONSTRAINT_MAT_H__
#define __FAUST_CONSTRAINT_MAT_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"
#include "faust_mat.h"


//template<typename parameter_type>
class faust_constraint_mat : public faust_constraint_generic
{
   public:
      faust_constraint_mat(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_, 
         const int nb_rows_, 
         const int nb_cols_);

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_, 
         const int nb_rows_, 
         const int nb_cols_,
         const faust_mat parameter_);

      faust_constraint_mat(const faust_constraint_mat& constraint_);

      faust_mat getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;
 
      ~faust_constraint_mat(){};

   private:
      // parameter of constraint
      faust_mat parameter;
    
};

#endif
