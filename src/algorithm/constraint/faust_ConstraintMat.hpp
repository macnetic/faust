/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __FAUST_CONSTRAINT_MAT_HPP__
#define __FAUST_CONSTRAINT_MAT_HPP__

//#include "faust_ConstraintMat.h"
#include <iostream>
#include <cstdlib>
#include "faust_MatDense.h"
#include "faust_exception.h"

template<typename FPP,Device DEVICE>
const char * Faust::ConstraintMat<FPP,DEVICE>::m_className="Faust::ConstraintMat<FPP,DEVICE>::";

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat() :
   Faust::ConstraintGeneric<FPP,DEVICE>()
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraintName_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraintName_,
         nbRows_,
         nbCols_)
{
   set_default_parameter();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraintName_,
   const Faust::MatDense<FPP,DEVICE> defaultParameter_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraintName_,
         nbRows_,
         nbCols_),
         m_parameter(defaultParameter_)
{
   check_constraint_name();
}

template<typename FPP,Device DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const Faust::ConstraintMat<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric<FPP,DEVICE>(
         constraint_.constraintName,
         constraint_.nbRows,
         constraint_.nbCols),
         m_parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
         break;
      case CONSTRAINT_NAME_SUPP:
         break;
      default:
         handleError(m_className," cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::set_default_parameter()
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
         m_parameter.setZeros();
         break;
      case CONSTRAINT_NAME_SUPP:
         m_parameter.setZeros();
         break;
      default:
         handleError(m_className,"set_default_parameter : cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}


template<typename FPP,Device DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
         mat=m_parameter;
         break;
      case CONSTRAINT_NAME_SUPP:
         Faust::prox_supp(mat,m_parameter);
         break;
      default:
         handleError(m_className,"project : invalid constraint_name");
         break;
   }
}



#endif
