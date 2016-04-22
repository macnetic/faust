#ifndef __FAUST_CONSTRAINT_INT_H__
#define __FAUST_CONSTRAINT_INT_H__

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
   template<typename FPP> class faust_cu_mat;
#else
   template<typename FPP> class faust_mat;
#endif


template<typename FPP>
class faust_constraint_int : public faust_constraint_generic<FPP>
{

   	#ifdef __COMPILE_GPU__
    		typedef faust_cu_mat<FPP> faust_matrix ;
	#else
    		typedef faust_mat<FPP> faust_matrix ;
	#endif		
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
      virtual void project(faust_matrix & mat)const;	
 
      ~faust_constraint_int(){};

   private:
      // parameter of constraint
      faust_unsigned_int parameter;
	  static const char * class_name;
    
};

#include "faust_constraint_int.hpp"

#endif
