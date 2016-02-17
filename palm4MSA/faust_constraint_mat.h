#ifndef __FAUST_CONSTRAINT_MAT_H__
#define __FAUST_CONSTRAINT_MAT_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"

#ifdef __COMPILE_GPU__
   #include "faust_cu_mat.h"
#else
   #include "faust_mat.h"
#endif


#ifdef __COMPILE_GPU__
   template<typename T> class faust_cu_mat;
#else
   template<typename T> class faust_mat;
#endif


//template<typename parameter_type>
template<typename T>
class faust_constraint_mat : public faust_constraint_generic
{
#ifdef __COMPILE_GPU__
   typedef faust_cu_mat<T> faust_matrix ;
#else
   typedef faust_mat<T> faust_matrix ;
#endif

   public:
      faust_constraint_mat(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_, 
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_mat(
         const faust_constraint_name& constraint_name_,  
         const faust_matrix parameter_,
         const faust_unsigned_int nb_rows_, 
         const faust_unsigned_int nb_cols_);

      faust_constraint_mat(const faust_constraint_mat& constraint_);

      faust_matrix getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;
 
      ~faust_constraint_mat(){};

   private:
      // parameter of constraint
      faust_matrix parameter;
	  static const char * class_name;
    
};

#include "faust_constraint_mat.hpp"

#endif
