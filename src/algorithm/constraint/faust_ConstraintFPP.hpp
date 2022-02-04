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
#ifndef __FAUST_CONSTRAINT_FPP_HPP__
#define __FAUST_CONSTRAINT_FPP_HPP__


#include <cstdlib>
#include "faust_exception.h"

template<typename FPP,FDevice DEVICE, typename FPP2>
const char * Faust::ConstraintFPP<FPP,DEVICE,FPP2>::m_className ="Faust::ConstraintFPP<FPP,DEVICE,FPP2>::";

template<typename FPP,FDevice DEVICE, typename FPP2>
Faust::ConstraintFPP<FPP,DEVICE,FPP2>::ConstraintFPP() :
   Faust::ConstraintGeneric()
{
   set_default_parameter();
}

template<typename FPP,FDevice DEVICE, typename FPP2>
Faust::ConstraintFPP<FPP,DEVICE,FPP2>::ConstraintFPP(
   const faust_constraint_name& constraintName_,
   const int nbRows_,
   const int nbCols_) :
      Faust::ConstraintGeneric(
         constraintName_,
         nbRows_,
         nbCols_)
{
   set_default_parameter();
}

template<typename FPP,FDevice DEVICE, typename FPP2>
Faust::ConstraintFPP<FPP,DEVICE,FPP2>::ConstraintFPP(
   const faust_constraint_name& constraintName_,
   const FPP2 defaultParameter_,
   const int nbRows_,
   const int nbCols_) :
      Faust::ConstraintGeneric(
         constraintName_,
         nbRows_,
         nbCols_),
         m_parameter(defaultParameter_)
{
   check_constraint_name();
}

template<typename FPP,FDevice DEVICE, typename FPP2>
Faust::ConstraintFPP<FPP,DEVICE,FPP2>::ConstraintFPP(
   const Faust::ConstraintFPP<FPP,DEVICE,FPP2>& constraint_) :
      Faust::ConstraintGeneric(
         constraint_.constraintName,
         constraint_.nbRows,
         constraint_.nbCols),
         m_parameter(constraint_.parameter)
{
   check_constraint_name();
}



template<typename FPP,FDevice DEVICE, typename FPP2>
void Faust::ConstraintFPP<FPP,DEVICE,FPP2>::check_constraint_name()const
{
   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_NORMCOL:
         break;
      case CONSTRAINT_NAME_NORMLIN:
         break;
      default:
         handleError(m_className,"check_constraint_name : cannot create Faust::ConstraintFPP objet from an faust_constraint object with this constraint_name");
         break;
   }
}

template<typename FPP,FDevice DEVICE, typename FPP2>
void Faust::ConstraintFPP<FPP,DEVICE,FPP2>::set_default_parameter()
{

   switch (this->m_constraintName)
   {
      case CONSTRAINT_NAME_NORMCOL:
         m_parameter = 0.0;
         break;
      case CONSTRAINT_NAME_NORMLIN:
         m_parameter = 0.0;
         break;
      default:
         handleError(m_className,"set_default_parameter : cannot create Faust::ConstraintFPP objet from an faust_constraint object with this->constraint_name");
         break;
   }
}

template<typename FPP,FDevice DEVICE,typename FPP2>
void Faust::ConstraintFPP<FPP,DEVICE,FPP2>::project(Faust::MatDense<FPP,DEVICE> & mat)const
{
	// no normalization on the whole matrix for this kind of constraint/prox (normalizing attribute ignored)
	switch (this->m_constraintName)
   	{
      		case CONSTRAINT_NAME_NORMCOL:
         		Faust::prox_normcol<FPP, FPP2>(mat,m_parameter);
         	break;
      		case CONSTRAINT_NAME_NORMLIN:
         		Faust::prox_normlin<FPP, FPP2>(mat,m_parameter);
         	break;
      		default:
         		handleError(m_className,"project : invalid constraint name");
         	break;
   }

}

template<typename FPP,FDevice DEVICE, typename FPP2>
Faust::MatGeneric<FPP,DEVICE>* Faust::ConstraintFPP<FPP,DEVICE, FPP2>::project_gen(Faust::MatDense<FPP,DEVICE> & mat) const
{
	switch (this->m_constraintName)
   	{
      		case CONSTRAINT_NAME_NORMCOL:
         		return Faust::prox_normcol_gen(mat, m_parameter);
         	break;
      		case CONSTRAINT_NAME_NORMLIN:
         		return Faust::prox_normlin_gen(mat, m_parameter);
         	break;
      		default:
         		handleError(m_className,"project : invalid constraint name");
         	break;
   }
}


template<typename FPP,FDevice DEVICE, typename FPP2>
const char* Faust::ConstraintFPP<FPP,DEVICE, FPP2>::get_type() const
{
	return "FAUST_REAL";
}

template<typename FPP,FDevice DEVICE, typename FPP2>
void Faust::ConstraintFPP<FPP,DEVICE,FPP2>::Display() const
{
	Faust::ConstraintGeneric::Display();
	std::cout<<" parameter :"<< get_parameter() << std::endl;
}

#endif

