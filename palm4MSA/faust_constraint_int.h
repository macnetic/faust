#ifndef __FAUST_CONSTRAINT_INT_H__
#define __FAUST_CONSTRAINT_INT_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"

//template<typename parameter_type>
class faust_constraint_int : public faust_constraint_generic
{
   public:
      faust_constraint_int(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_int(
         const faust_constraint_name& constraint_name_, 
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_int(
         const faust_constraint_name& constraint_name_, 
         const faust_unsigned_int parameter_,
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_int(const faust_constraint_int& constraint_);

      faust_unsigned_int getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;

 
      ~faust_constraint_int(){};

   private:
      // parameter of constraint
      faust_unsigned_int parameter;
    
};

#endif
