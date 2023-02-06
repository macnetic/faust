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
#ifndef __FAUST_CONSTRAINT_MAT_HPP__
#define __FAUST_CONSTRAINT_MAT_HPP__

//#include "faust_ConstraintMat.h"
#include <iostream>
#include <cstdlib>
#include "faust_MatDense.h"
#include "faust_exception.h"

template<typename FPP,FDevice DEVICE>
const char * Faust::ConstraintMat<FPP,DEVICE>::m_className="Faust::ConstraintMat<FPP,DEVICE>::";

template<typename FPP,FDevice DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat() :
   Faust::ConstraintGeneric()
{
   set_default_parameter();
}

template<typename FPP,FDevice DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraintName_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric(
         constraintName_,
         nbRows_,
         nbCols_)
{
   set_default_parameter();
}

template<typename FPP,FDevice DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const faust_constraint_name& constraintName_,
   const Faust::MatDense<FPP,DEVICE> defaultParameter_,
   const faust_unsigned_int nbRows_,
   const faust_unsigned_int nbCols_) :
      Faust::ConstraintGeneric(
         constraintName_,
         nbRows_,
         nbCols_),
         m_parameter(defaultParameter_)
{
   check_constraint_name();
}

template<typename FPP,FDevice DEVICE>
Faust::ConstraintMat<FPP,DEVICE>::ConstraintMat(
   const Faust::ConstraintMat<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric(
         constraint_.constraintName,
         constraint_.nbRows,
         constraint_.nbCols),
         m_parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename FPP,FDevice DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
      case CONSTRAINT_NAME_SUPP:
	  case CONSTRAINT_NAME_CIRC:
	  case CONSTRAINT_NAME_ANTICIRC:
	  case CONSTRAINT_NAME_TOEPLITZ:
	  case CONSTRAINT_NAME_HANKEL:
	  case CONSTRAINT_NAME_BLKDIAG:
	  case CONSTRAINT_NAME_ID:
			 break;
      default:
         handleError(m_className," cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename FPP,FDevice DEVICE>
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
	  case CONSTRAINT_NAME_ID:
		 // nothing to do
		 break;
      default:
         handleError(m_className,"set_default_parameter : cannot create Faust::ConstraintMat objet from an faust_constraint object with this constraint_name");
         break;
   }
}


template<typename FPP,FDevice DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
         //mat=m_parameter;
		 Faust::prox_const(mat, m_parameter, normalizing, pos);
         break;
      case CONSTRAINT_NAME_SUPP:
         Faust::prox_supp(mat,m_parameter, normalizing, pos);
         break;
	  case CONSTRAINT_NAME_TOEPLITZ:
		 Faust::prox_toeplitz(mat, normalizing, pos);
		 break;
	  case CONSTRAINT_NAME_CIRC:
		 Faust::prox_circ(mat, normalizing, pos);
		 break;
	  case CONSTRAINT_NAME_ANTICIRC:
		 Faust::prox_anticirc(mat, normalizing, pos);
		 break;
	  case CONSTRAINT_NAME_HANKEL:
		 Faust::prox_hankel(mat, normalizing, pos);
		 break;
	  case CONSTRAINT_NAME_BLKDIAG:
		 Faust::prox_blockdiag(mat, m_parameter, normalizing, pos);
		 break;
	  case CONSTRAINT_NAME_ID:
		 Faust::prox_id(mat, normalizing, pos);
		 break;
      default:
         handleError(m_className,"project : invalid constraint_name");
         break;
   }
}

template<typename FPP, FDevice DEVICE>
Faust::MatGeneric<FPP,DEVICE>* Faust::ConstraintMat<FPP,DEVICE>::project_gen(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_CONST:
         //mat=m_parameter;
		return  Faust::prox_const_gen(mat, m_parameter, normalizing, pos);
      case CONSTRAINT_NAME_ID:
		return  Faust::prox_id_gen(mat, normalizing, pos);
      case CONSTRAINT_NAME_SUPP:
         return Faust::prox_supp_gen(mat,m_parameter, normalizing, pos);
	  case CONSTRAINT_NAME_TOEPLITZ:
		 return Faust::prox_toeplitz_gen(mat, normalizing, pos);
	  case CONSTRAINT_NAME_CIRC:
		 return Faust::prox_circ_gen(mat, normalizing, pos);
	  case CONSTRAINT_NAME_ANTICIRC:
		 return Faust::prox_anticirc_gen(mat, normalizing, pos);
	  case CONSTRAINT_NAME_HANKEL:
		 return Faust::prox_hankel_gen(mat, normalizing, pos);
	  case CONSTRAINT_NAME_BLKDIAG:
		 return Faust::prox_blockdiag_gen(mat, m_parameter, normalizing, pos);
      default:
         handleError(m_className,"project : invalid constraint_name");
         break;
   }
}
template<typename FPP,FDevice DEVICE>
void Faust::ConstraintMat<FPP,DEVICE>::Display() const
{
	Faust::ConstraintGeneric::Display();
	std::cout<<" parameter :";
	get_parameter().Display();
}

template<typename FPP,FDevice DEVICE>
const char* Faust::ConstraintMat<FPP,DEVICE>::get_type() const
{
	return "FAUST_MAT";
}

#endif

