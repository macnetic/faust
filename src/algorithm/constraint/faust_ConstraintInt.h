#ifndef __FAUST_CONSTRAINT_INT_H__
#define __FAUST_CONSTRAINT_INT_H__

#include "faust_ConstraintGeneric.h"
#include "faust_constant.h"



#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
   #include "faust_prox_gpu.h"
#else
   #include "faust_MatDense.h"
   #include "faust_prox.h"
#endif


template<typename FPP,Device DEVICE> class MatDense;

//! \class ConstraintInt
//! \brief Contains the integer constraint parameters for the hierarchical factorization. <br>
//! See the parent class Faust::ConstraintGeneric for more detail. <br>


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template<typename FPP,Device DEVICE> class ConstraintGeneric;

    template<typename FPP,Device DEVICE>
    class ConstraintInt : public Faust::ConstraintGeneric<FPP,DEVICE>
    {


       public:
          ConstraintInt(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

          ConstraintInt(
             const faust_constraint_name& constraintName_,
             const faust_unsigned_int nbRows_,
             const faust_unsigned_int nbCols_);

          ConstraintInt(
             const faust_constraint_name& constraintName_,
             const faust_unsigned_int parameter_,
             const faust_unsigned_int nbRows_,
             const faust_unsigned_int nbCols_);

          ConstraintInt(const ConstraintInt& constraint_);

          faust_unsigned_int get_parameter() const {return m_parameter;};

          virtual void set_default_parameter();
          virtual void check_constraint_name()const;
          virtual void project(Faust::MatDense<FPP,DEVICE> & mat)const;

          ~ConstraintInt(){};

       private:
          // parameter of constraint
          faust_unsigned_int m_parameter;
          static const char * m_className;

    };
}

#include "faust_ConstraintInt.hpp"

#endif
