#ifndef __FAUST_CONSTRAINT_MAT_H__
#define __FAUST_CONSTRAINT_MAT_H__

#include "faust_ConstraintGeneric.h"
#include "faust_constant.h"

#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
   #include "faust_prox_gpu.h"
#else
   #include "faust_MatDense.h"
   #include "faust_prox.h"
#endif
//! \class Faust::ConstraintMat
//!    \brief Contains the matrix dense constraint parameters for the hierarchical factorization. <br>
//!    See the parent class Faust::ConstraintGeneric for more detail. <br>


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template<typename FPP,Device DEVICE> class ConstraintGeneric;
    template<typename FPP,Device DEVICE> class MatDense;

    template<typename FPP,Device DEVICE>
    class ConstraintMat : public Faust::ConstraintGeneric<FPP,DEVICE>
    {

       public:
          ConstraintMat(); // ajouter parametre de contrainte par defaut (a voir avec Luc)

          ConstraintMat(
             const faust_constraint_name& constraintName_,
             const faust_unsigned_int nbRows_,
             const faust_unsigned_int nbCols_);

          ConstraintMat(
             const faust_constraint_name& constraintName_,
             const Faust::MatDense<FPP,DEVICE> parameter_,
             const faust_unsigned_int nbRows_,
             const faust_unsigned_int nbCols_);

          ConstraintMat(const ConstraintMat& constraint_);

          Faust::MatDense<FPP,DEVICE> get_parameter() const {return m_parameter;};

          virtual void set_default_parameter();
          virtual void check_constraint_name()const;
          virtual void project(Faust::MatDense<FPP,DEVICE> & mat)const;
          ~ConstraintMat(){};

       private:
          // parameter of constraint
          Faust::MatDense<FPP,DEVICE> m_parameter;
          static const char * m_className;

    };
}

#include "faust_ConstraintMat.hpp"

#endif
