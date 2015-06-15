#ifndef __FAUST_CONSTRAINT_REAL_H__
#define __FAUST_CONSTRAINT_REAL_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"

//template<typename parameter_type>
class faust_constraint_real : public faust_constraint_generic
{
   public:
      faust_constraint_real(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_real(
         const faust_constraint_name& constraint_name_, 
         const int nb_rows_, 
         const int nb_cols_);

      faust_constraint_real(
         const faust_constraint_name& constraint_name_,  
         const faust_real parameter_,
         const int nb_rows_, 
         const int nb_cols_);

      faust_constraint_real(const faust_constraint_real& constraint_);

      faust_real getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;
 
      ~faust_constraint_real(){};

   private:
      // parameter of constraint
      faust_real parameter;


      
};

#endif
