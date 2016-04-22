#ifndef __FAUST_CONSTRAINT_REAL_H__
#define __FAUST_CONSTRAINT_REAL_H__

#include "faust_constraint_generic.h"
#include "faust_constant.h"


#ifdef __COMPILE_GPU__
   #include "faust_cu_mat.h"
   #include "prox_cu.h"
#else
   #include "faust_mat.h"
   #include "prox.h"
#endif

#ifdef __COMPILE_GPU__
   template<typename T> class faust_cu_mat;
#else
   template<typename T> class faust_mat;
#endif


template<typename T> 
class faust_constraint_real : public faust_constraint_generic<T>
{

	#ifdef __COMPILE_GPU__
    		typedef faust_cu_mat<T> faust_matrix ;
	#else
    		typedef faust_mat<T> faust_matrix ;
	#endif	
	
   public:
      faust_constraint_real(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

      faust_constraint_real(
         const faust_constraint_name& constraint_name_, 
         const int nb_rows_, 
         const int nb_cols_);

      faust_constraint_real(
         const faust_constraint_name& constraint_name_,  
         const T parameter_,
         const int nb_rows_, 
         const int nb_cols_);

      faust_constraint_real(const faust_constraint_real& constraint_);

      T getParameter() const {return parameter;};
      
      virtual void set_default_parameter();
      virtual void check_constraint_name()const;
      virtual void project(faust_matrix & mat)const;
      ~faust_constraint_real(){};

   private:
      // parameter of constraint
      T parameter;
	  static const char * class_name;


      
};

#include "faust_constraint_real.hpp"

#endif
