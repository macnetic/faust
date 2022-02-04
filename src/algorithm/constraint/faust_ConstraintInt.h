/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_CONSTRAINT_INT_H__
#define __FAUST_CONSTRAINT_INT_H__

#include "faust_ConstraintGeneric.h"
#include "faust_constant.h"



#include "faust_MatDense.h"
#include "faust_prox.h"
#if USE_GPU_MOD
#include "faust_prox_gpu.h"
#include "faust_MatDense_gpu.h"
#endif
#include "faust_prox_gen.h"


template<typename FPP,FDevice DEVICE> class MatDense;

//! \class ConstraintInt
//! \brief Contains the integer constraint parameters for the hierarchical factorization. <br>
//! See the parent class Faust::ConstraintGeneric for more detail. <br>


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    class ConstraintGeneric;

    template<typename FPP,FDevice DEVICE>
    class ConstraintInt : public ConstraintGeneric
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
		  const char* get_type() const;
          virtual void project(Faust::MatDense<FPP,DEVICE> & mat)const;
          virtual MatGeneric<FPP,DEVICE>* project_gen(Faust::MatDense<FPP,DEVICE> & mat) const;
		  virtual void Display() const;
          ~ConstraintInt(){};

       private:
          // parameter of constraint
          faust_unsigned_int m_parameter;
          static const char * m_className;

    };
}

#include "faust_ConstraintInt.hpp"

#endif
