/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
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
#ifndef __FAUST_CONSTRAINT_GENERIC_HPP__
#define __FAUST_CONSTRAINT_GENERIC_HPP__

#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include <typeinfo>
#include "faust_exception.h"

#include "faust_ConstraintType.h"

//modif AL AL
//template<typename FPP,Device DEVICE> class ConstraintInt;

template<typename FPP,Device DEVICE> class ConstraintFPP;
template<typename FPP,Device DEVICE> class ConstraintMat;
template<typename FPP,Device DEVICE> class faust_ConstraintType;

template<typename FPP,Device DEVICE>
const char * Faust::ConstraintGeneric<FPP,DEVICE>::m_className="Faust::ConstraintGeneric::";

template<typename FPP,Device DEVICE>
const faust_constraint_name Faust::ConstraintGeneric<FPP,DEVICE>::get_constraint_type() const
{
   return m_constraintName;
}





template<typename FPP,Device DEVICE>
const char* Faust::ConstraintGeneric<FPP,DEVICE>::get_constraint_name()const
{
   switch(m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         return "CONSTRAINT_NAME_SP";
      case CONSTRAINT_NAME_SPCOL:
         return "CONSTRAINT_NAME_SPCOL";
      case CONSTRAINT_NAME_SPLIN:
         return "CONSTRAINT_NAME_SPLIN";
      case CONSTRAINT_NAME_NORMCOL:
         return "CONSTRAINT_NAME_NORMCOL";
      case CONSTRAINT_NAME_SPLINCOL:
         return "CONSTRAINT_NAME_SPLINCOL";
      case CONSTRAINT_NAME_CONST:
         return "CONSTRAINT_NAME_CONST";
      case CONSTRAINT_NAME_SP_POS:
         return "CONSTRAINT_NAME_SP_POS";
      case CONSTRAINT_NAME_BLKDIAG:
         return "CONSTRAINT_NAME_BLKDIAG";
      case CONSTRAINT_NAME_SUPP:
         return "CONSTRAINT_NAME_SUPP";
      case CONSTRAINT_NAME_NORMLIN:
         return "CONSTRAINT_NAME_NORMLIN";
      default:
         return "unknown constraint name";
   }
}


template<typename FPP,Device DEVICE>
const char*  Faust::ConstraintGeneric<FPP,DEVICE>::get_type() const
{
   switch(m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
           handleError(m_className,"get_type : unknown type parameter");
		   }
      case CONSTRAINT_NAME_SPLINCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
		 }
      case CONSTRAINT_NAME_CONST:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SP_POS:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_BLKDIAG:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SUPP:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
			handleError(m_className,"get_type : unknown type parameter");
		   }
      default:
			handleError(m_className,"get_type : unknown constraint type ");
   }
}


template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::is_constraint_parameter_int()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(m_className,"is_constraint_parameter_int : Unknown type of constraint");
		break;
	}
    return false;
}

template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::is_constraint_parameter_real()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(m_className,"is_constraint_parameter_real : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::is_constraint_parameter_mat()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSplincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeConst)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSpPos)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeSupp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::ConstraintTypeNormlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(m_className,"is_constraint_parameter_mat : Unknown type of constraint");
		break;
	}
        return false;
}


#endif
