#include "faust_constraint_generic.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_type.h"
#include <typeinfo>
#include "faust_exception.h"


const char * faust_constraint_generic::class_name="faust_constraint_generic::"; 

const faust_constraint_name faust_constraint_generic::getConstraintType() const 
{
   return constraint_name;
}

const char*  faust_constraint_generic::getType() const
{	
   switch(constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         if(typeid(constraint_type_sp)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_sp)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_sp)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPCOL:
         if(typeid(constraint_type_spcol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_spcol)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_spcol)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN:
         if(typeid(constraint_type_splin)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_splin)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_splin)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMCOL:
         if(typeid(constraint_type_normcol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_normcol)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_normcol)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
           handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_SPLINCOL:
         if(typeid(constraint_type_splincol)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_splincol)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_splincol)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_L0PEN:
         if(typeid(constraint_type_l0pen)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_l0pen)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_l0pen)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_L1PEN:
         if(typeid(constraint_type_l1pen)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_l1pen)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_l1pen)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
           handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_CONST:
         if(typeid(constraint_type_const)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_const)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_const)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
            handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_WAV:
         if(typeid(constraint_type_wav)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_wav)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_wav)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SP_POS:
         if(typeid(constraint_type_sp_pos)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_sp_pos)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_sp_pos)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_BLKDIAG:
         if(typeid(constraint_type_blkdiag)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_blkdiag)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_blkdiag)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SPLIN_TEST:
         if(typeid(constraint_type_splin_test)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_splin_test)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_splin_test)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_SUPP:
         if(typeid(constraint_type_supp)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_supp)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_supp)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
				handleError(class_name,"getType : unknown type parameter");
			}
      case CONSTRAINT_NAME_NORMLIN:
         if(typeid(constraint_type_normlin)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_normlin)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_normlin)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
			handleError(class_name,"getType : unknown type parameter");
		   }
      case CONSTRAINT_NAME_TOEPLITZ:
         if(typeid(constraint_type_toeplitz)==typeid(faust_constraint_int))
            return "INT";
         else if(typeid(constraint_type_toeplitz)==typeid(faust_constraint_real))
            return "FAUST_REAL";
         else if(typeid(constraint_type_toeplitz)==typeid(faust_constraint_mat))
            return "FAUST_MAT";
         else{
			handleError(class_name,"getType : unknown type parameter");
		   }
      default:
			handleError(class_name,"getType : unknown constraint type ");
   }
}

bool faust_constraint_generic::isConstraintParameterInt()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(constraint_type_sp)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(constraint_type_spcol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(constraint_type_splin)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(constraint_type_normcol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(constraint_type_splincol)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(constraint_type_l0pen)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(constraint_type_l1pen)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(constraint_type_const)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(constraint_type_wav)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(constraint_type_sp_pos)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(constraint_type_blkdiag)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(constraint_type_splin_test)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(constraint_type_supp)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(constraint_type_normlin)==typeid(faust_constraint_int)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(constraint_type_toeplitz)==typeid(faust_constraint_int)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterInt : Unknown type of constraint");
		break;
	}
        return false;
}

bool faust_constraint_generic::isConstraintParameterReal()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(constraint_type_sp)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(constraint_type_spcol)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(constraint_type_splin)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(constraint_type_normcol)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(constraint_type_splincol)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(constraint_type_l0pen)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(constraint_type_l1pen)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(constraint_type_const)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(constraint_type_wav)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(constraint_type_sp_pos)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(constraint_type_blkdiag)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(constraint_type_splin_test)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(constraint_type_supp)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(constraint_type_normlin)==typeid(faust_constraint_real)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(constraint_type_toeplitz)==typeid(faust_constraint_real)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterReal : Unknown type of constraint");
		break;
	}
        return false;
}

bool faust_constraint_generic::isConstraintParameterMat()const
{
	switch(constraint_name)
	{
		case CONSTRAINT_NAME_SP:
			return (typeid(constraint_type_sp)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SPCOL:
			return (typeid(constraint_type_spcol)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN:
			return (typeid(constraint_type_splin)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_NORMCOL:
			return (typeid(constraint_type_normcol)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SPLINCOL:
			return (typeid(constraint_type_splincol)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_L0PEN:
			return (typeid(constraint_type_l0pen)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_L1PEN:
			return (typeid(constraint_type_l1pen)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_CONST:
			return (typeid(constraint_type_const)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_WAV:
			return (typeid(constraint_type_wav)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SP_POS:
			return (typeid(constraint_type_sp_pos)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_BLKDIAG:
			return (typeid(constraint_type_blkdiag)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SPLIN_TEST:
			return (typeid(constraint_type_splin_test)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_SUPP:
			return (typeid(constraint_type_supp)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_NORMLIN:
			return (typeid(constraint_type_normlin)==typeid(faust_constraint_mat)?true:false);
		break;
		case CONSTRAINT_NAME_TOEPLITZ:
			return (typeid(constraint_type_toeplitz)==typeid(faust_constraint_mat)?true:false);
		break;
		default:
			handleError(class_name,"isConstraintParameterMat : Unknown type of constraint");
		break;
	}
        return false;
}




const char* faust_constraint_generic::get_constraint_name()const
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
      case CONSTRAINT_NAME_L0PEN:
         return "CONSTRAINT_NAME_L0PEN";
      case CONSTRAINT_NAME_L1PEN:
         return "CONSTRAINT_NAME_L1PEN";
      case CONSTRAINT_NAME_CONST:
         return "CONSTRAINT_NAME_CONST";
      case CONSTRAINT_NAME_WAV:
         return "CONSTRAINT_NAME_WAV";
      case CONSTRAINT_NAME_SP_POS:
         return "CONSTRAINT_NAME_SP_POS";
      case CONSTRAINT_NAME_BLKDIAG:
         return "CONSTRAINT_NAME_BLKDIAG";
      case CONSTRAINT_NAME_SPLIN_TEST:
         return "CONSTRAINT_NAME_SPLIN_TEST";
      case CONSTRAINT_NAME_SUPP:
         return "CONSTRAINT_NAME_SUPP";
      case CONSTRAINT_NAME_NORMLIN:
         return "CONSTRAINT_NAME_NORMLIN";
      case CONSTRAINT_NAME_TOEPLITZ:
         return "CONSTRAINT_NAME_TOEPLITZ";
      default:
         return "unknown constraint name";
   }
}


