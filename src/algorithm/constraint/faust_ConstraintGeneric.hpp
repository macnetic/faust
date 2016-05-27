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
const char * Faust::ConstraintGeneric<FPP,DEVICE>::class_name="Faust::ConstraintGeneric::";

template<typename FPP,Device DEVICE>
const faust_constraint_name Faust::ConstraintGeneric<FPP,DEVICE>::getConstraintType() const
{
   return constraint_name;
}





template<typename FPP,Device DEVICE>
const char* Faust::ConstraintGeneric<FPP,DEVICE>::get_constraint_name()const
{
   switch(constraint_name)
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
const char*  Faust::ConstraintGeneric<FPP,DEVICE>::getType() const
{
   switch(constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
           handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_SPLINCOL:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
		 }
      case CONSTRAINT_NAME_CONST:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SP_POS:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_BLKDIAG:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SUPP:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMLIN:
         if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>))
            return "INT";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>))
            return "FAUST_REAL";
         else if(typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>))
            return "FAUST_MAT";
         else{
			handleError(class_name,"getType : unknown type parameter");
		   }
      default:
			handleError(class_name,"getType : unknown constraint type ");
   }
}


template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::isConstraintParameterInt()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintInt<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterInt : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::isConstraintParameterReal()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintFPP<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterReal : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename FPP,Device DEVICE>
bool Faust::ConstraintGeneric<FPP,DEVICE>::isConstraintParameterMat()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_spcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normcol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_splincol)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_const)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_sp_pos)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_blkdiag)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_supp)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  ConstraintType<FPP,DEVICE>::constraint_type_normlin)==typeid(Faust::ConstraintMat<FPP,DEVICE>)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterMat : Unknown type of constraint");
		break;
	}
        return false;
}


#endif
