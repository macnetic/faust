#ifndef __FAUST_CONSTRAINT_FPP_H__
#define __FAUST_CONSTRAINT_FPP_H__

#include "faust_ConstraintGeneric.h"
#include "faust_constant.h"


#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
   #include "prox_gpu.h"
#else
   #include "faust_MatDense.h"
   #include "prox.h"
#endif



/*!
 \class Faust::ConstraintFPP
 \brief Contains the Floating Point Precision constraint parameters defined by template FPP for the hierarchical factorization. <br>
FPP can be float or double precision.
See the parent class Faust::ConstraintGeneric for more detail. <br>
*/

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template<typename FPP,Device DEVICE> class ConstraintGeneric;
    template<typename FPP,Device DEVICE> class MatDense;

    template<typename FPP,Device DEVICE>
    class ConstraintFPP : public Faust::ConstraintGeneric<FPP,DEVICE>
    {



       public:
          ConstraintFPP(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

          ConstraintFPP(
             const faust_constraint_name& constraint_name_,
             const int nb_rows_,
             const int nb_cols_);

          ConstraintFPP(
             const faust_constraint_name& constraint_name_,
             const FPP parameter_,
             const int nb_rows_,
             const int nb_cols_);

          ConstraintFPP(const ConstraintFPP& constraint_);

          FPP getParameter() const {return parameter;};

          virtual void set_default_parameter();
          virtual void check_constraint_name()const;
          virtual void project(Faust::MatDense<FPP,DEVICE> & mat)const;
          ~ConstraintFPP(){};

       private:
          // parameter of constraint
          FPP parameter;
          static const char * class_name;



    };
}

#include "faust_ConstraintFPP.hpp"

#endif
