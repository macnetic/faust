#include "faust_constraint_generic.h"


std::string getConstraintType(faust_constraint_name cons_name)  
{
   switch(cons_name)
   {
      case CONSTRAINT_NAME_SP:
         return "INT";
      case CONSTRAINT_NAME_SPCOL:
         return "INT";
      case CONSTRAINT_NAME_SPLIN:
         return "INT";
      case CONSTRAINT_NAME_NORMCOL:
         return "FAUST_REAL";
      case CONSTRAINT_NAME_SPLINCOL:
         return "INT";
      case CONSTRAINT_NAME_L0PEN:
         return "INT";
      case CONSTRAINT_NAME_L1PEN:
         return "INT";
      case CONSTRAINT_NAME_CONST:
         return "FAUST_MAT";
      case CONSTRAINT_NAME_WAV:
         return "INT";
      case CONSTRAINT_NAME_SP_POS:
         return "INT";
      case CONSTRAINT_NAME_BLKDIAG:
         return "INT";
      case CONSTRAINT_NAME_SPLIN_TEST:
         return "INT";
      case CONSTRAINT_NAME_SUPP:
         return "FAUST_MAT";
      case CONSTRAINT_NAME_NORMLIN:
         return "FAUST_REAL";
      case CONSTRAINT_NAME_TOEPLITZ:
         return "INT";
      default:
         return "unknown constraint type";
   }
}

char*  faust_constraint_generic::getType() const
{	
	
   switch(constraint_name)
   {
      case CONSTRAINT_NAME_SP:
         return "INT";
      case CONSTRAINT_NAME_SPCOL:
         return "INT";
      case CONSTRAINT_NAME_SPLIN:
         return "INT";
      case CONSTRAINT_NAME_NORMCOL:
         return "FAUST_REAL";
      case CONSTRAINT_NAME_SPLINCOL:
         return "INT";
      case CONSTRAINT_NAME_L0PEN:
         return "INT";
      case CONSTRAINT_NAME_L1PEN:
         return "INT";
      case CONSTRAINT_NAME_CONST:
         return "FAUST_MAT";
      case CONSTRAINT_NAME_WAV:
         return "INT";
      case CONSTRAINT_NAME_SP_POS:
         return "INT";
      case CONSTRAINT_NAME_BLKDIAG:
         return "INT";
      case CONSTRAINT_NAME_SPLIN_TEST:
         return "INT";
      case CONSTRAINT_NAME_SUPP:
         return "FAUST_MAT";
      case CONSTRAINT_NAME_NORMLIN:
         return "FAUST_REAL";
      case CONSTRAINT_NAME_TOEPLITZ:
         return "INT";
      default:
         return "unknown constraint type";
   }
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


