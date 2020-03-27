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

#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include <typeinfo>
#include "faust_exception.h"

#include "faust_ConstraintType.h"

//modif AL AL
//template<typename FPP,FDevice DEVICE> class ConstraintInt;

template<typename FPP,FDevice DEVICE,typename FPP2> class ConstraintFPP;
template<typename FPP,FDevice DEVICE> class ConstraintMat;
template<typename FPP,FDevice DEVICE,typename FPP2> class faust_ConstraintType;


template<typename FPP, FDevice DEVICE, typename FPP2>
const char*  Faust::ConstraintGeneric::get_type() const
{
   switch(m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
           handleError(m_className,"get_type : unknown type parameter");
		   }
      case CONSTRAINT_NAME_SPLINCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
		 }
      case CONSTRAINT_NAME_CONST:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SP_POS:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_BLKDIAG:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_SUPP:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(m_className,"get_type : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
			handleError(m_className,"get_type : unknown type parameter");
		   }
	  case CONSTRAINT_NAME_TOEPLITZ:
	  case CONSTRAINT_NAME_CIRC:
	  case CONSTRAINT_NAME_HANKEL:
		 return "FAUST_MAT";
      default:
			handleError(m_className,"get_type : unknown constraint type ");
   }
}


template<typename FPP,FDevice DEVICE, typename FPP2>
bool Faust::ConstraintGeneric::is_constraint_parameter_int()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
		case CONSTRAINT_NAME_HANKEL:
		case CONSTRAINT_NAME_CIRC:
			return false;
		default:
			handleError(m_className,"is_constraint_parameter_int : Unknown type of constraint");
		break;
	}
    return false;
}

template<typename FPP,FDevice DEVICE, typename FPP2>
bool Faust::ConstraintGeneric::is_constraint_parameter_real()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE,FPP2>)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
		case CONSTRAINT_NAME_HANKEL:
		case CONSTRAINT_NAME_CIRC:
			return false;
		default:
			handleError(m_className,"is_constraint_parameter_real : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename FPP,FDevice DEVICE, typename FPP2>
bool Faust::ConstraintGeneric::is_constraint_parameter_mat()const
{
	switch(m_constraintName)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSplincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeConst)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSpPos)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeBlkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeSupp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE,FPP2>::ConstraintTypeNormlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
		case CONSTRAINT_NAME_HANKEL:
		case CONSTRAINT_NAME_CIRC:
			return true;
		default:
			handleError(m_className,"is_constraint_parameter_mat : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename FPP,FDevice DEVICE, typename FPP2>
void Faust::ConstraintGeneric::project(Faust::MatDense<FPP, DEVICE>& mat) const {
	//unfortunately it's not possible to do template with virtual (pure or not) function
	// (it needs to be a template class with virtual function using template types to be possible)
	// (but we don't want that because it implies to pass all the templates types when instantiating the child classes
	// and this would be non-efficient because for example FPP2 is not needed by ConstraintMat or ConstraintInt)
	// (and we need FPP2 because the field to define the dense matrix has not be the same than the real type used to define
	// the real constraints, by the way FPP could even be a complex so we do really need  FPP2 even if it's heavy impl. Neverthless, there is a default template value type of float to lighten the use.)
//	std::cout << "Faust::ConstraintGeneric::project(mat) typeid="<< typeid(this).name() << std::endl;
	if(this->is_constraint_parameter_mat<FPP,DEVICE,FPP2>())
		dynamic_cast<const Faust::ConstraintMat<FPP,DEVICE>*>(this)->project(mat);
	else if(this->is_constraint_parameter_real<FPP,DEVICE,FPP2>()){
		dynamic_cast<const Faust::ConstraintFPP<FPP,DEVICE,FPP2>*>(this)->project(mat);
	}
	else if(this->is_constraint_parameter_int<FPP,DEVICE,FPP2>())
		dynamic_cast<const Faust::ConstraintInt<FPP,DEVICE>*>(this)->project(mat);
}


