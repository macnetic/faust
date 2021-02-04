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
#include "faust_prox.h"
#ifdef USE_GPU_MOD
#include "faust_prox_gpu.h"
#endif
#include "faust_ConstraintInt.h"
#include <iostream>
#include <cstdlib>
#include "faust_exception.h"


template<typename FPP,FDevice DEVICE>
const char * Faust::ConstraintInt<FPP,DEVICE>::m_className="Faust::ConstraintInt";

template<typename FPP,FDevice DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt() :
   Faust::ConstraintGeneric()
{
   set_default_parameter();
}

template<typename FPP,FDevice DEVICE>
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
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
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const faust_constraint_name& constraintName_,
   const faust_unsigned_int defaultParameter_,
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
Faust::ConstraintInt<FPP,DEVICE>::ConstraintInt(
   const Faust::ConstraintInt<FPP,DEVICE>& constraint_) :
      Faust::ConstraintGeneric(
         constraint_.constraintName,
         constraint_.nbRows,
         constraint_.nbCols),
         m_parameter(constraint_.parameter)
{
   check_constraint_name();
}


template<typename FPP,FDevice DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::check_constraint_name()const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
      case CONSTRAINT_NAME_SPCOL:
      case CONSTRAINT_NAME_SPLIN:
      case CONSTRAINT_NAME_SPLINCOL:
      case CONSTRAINT_NAME_SP_POS:
	  case CONSTRAINT_NAME_BLKDIAG: //TODO: shouldn't it be a matrix constraint?
	  case CONSTRAINT_NAME_SKPERM:
         break;
      default:
		handleError(m_className," cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with constraint_name");
         break;
   }
}

template<typename FPP,FDevice DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::set_default_parameter()
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPCOL:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLIN:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SP_POS:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         m_parameter = 0;
         break;
      case CONSTRAINT_NAME_SKPERM:
         m_parameter = 0;
         break;
      default:
		handleError(m_className,"set_default_parameter : cannot create Faust::ConstraintInt objet from an faust_constraint object with constraint with this constraint_name");
         break;
   }
}

template<typename FPP,FDevice DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::project(Faust::MatDense<FPP,DEVICE> & mat) const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         Faust::prox_sp(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPCOL:
         Faust::prox_spcol(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPLIN:
         Faust::prox_splin(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         Faust::prox_splincol(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SP_POS:
         Faust::prox_sp_pos(mat,m_parameter);
         break;
      case CONSTRAINT_NAME_SKPERM:
         Faust::prox_skperm(mat,m_parameter);
         break;
      default:
		handleError(m_className,"project : cannot project with this constraint name");
         break;
   }
}

template<typename FPP,FDevice DEVICE>
void Faust::ConstraintInt<FPP,DEVICE>::Display() const
{
	Faust::ConstraintGeneric::Display();
	std::cout<<" parameter :"<< get_parameter() << std::endl;
}


