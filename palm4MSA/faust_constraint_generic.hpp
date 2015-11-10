#include "faust_constraint_generic.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_type.h"
#include <typeinfo>
#include "faust_exception.h"

template<typename T> class faust_constraint_real;
template<typename T> class faust_constraint_mat;





template<typename T>
const char*  faust_constraint_generic::getType() const
{	
   switch(constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         if(typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPCOL:
         if(typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN:
         if(typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMCOL:
         if(typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
           handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_SPLINCOL:
         if(typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_L0PEN:
         if(typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_L1PEN:
         if(typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
           handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_CONST:
         if(typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_WAV:
         if(typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SP_POS:
         if(typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_BLKDIAG:
         if(typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN_TEST:
         if(typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SUPP:
         if(typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMLIN:
         if(typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
			handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_TOEPLITZ:
         if(typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_real<T>))
            return "FAUST_REAL";
         else if(typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_mat<T>))
            return "FAUST_MAT";
         else{
			handleError(class_name,"getType : unknown type parameter");
		   }
      default:
			handleError(class_name,"getType : unknown constraint type ");
   }
}


template<typename T>
bool faust_constraint_generic::isConstraintParameterInt()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_int)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterInt : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename T>
bool faust_constraint_generic::isConstraintParameterReal()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_real<T>)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_real<T>)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterReal : Unknown type of constraint");
		break;
	}
        return false;
}

template<typename T>
bool faust_constraint_generic::isConstraintParameterMat()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(typename  constraint_type<T>::constraint_type_sp)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_spcol)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_splin)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_normcol)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(typename  constraint_type<T>::constraint_type_splincol)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l0pen)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(typename  constraint_type<T>::constraint_type_l1pen)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(typename  constraint_type<T>::constraint_type_const)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(typename  constraint_type<T>::constraint_type_wav)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(typename  constraint_type<T>::constraint_type_sp_pos)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(typename  constraint_type<T>::constraint_type_blkdiag)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(typename  constraint_type<T>::constraint_type_splin_test)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(typename  constraint_type<T>::constraint_type_supp)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(typename  constraint_type<T>::constraint_type_normlin)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(typename  constraint_type<T>::constraint_type_toeplitz)==typeid(faust_constraint_mat<T>)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterMat : Unknown type of constraint");
		break;
	}
        return false;
}



