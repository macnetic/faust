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

template char const* faust_constraint_generic::getType<double>() const;
template bool faust_constraint_generic::isConstraintParameterInt<double>() const;
template bool faust_constraint_generic::isConstraintParameterReal<double>() const;
template bool faust_constraint_generic::isConstraintParameterMat<double>() const;

template char const* faust_constraint_generic::getType<float>() const;
template bool faust_constraint_generic::isConstraintParameterInt<float>() const;
template bool faust_constraint_generic::isConstraintParameterReal<float>() const;
template bool faust_constraint_generic::isConstraintParameterMat<float>() const;

bool isConstraintNameInt(const char * type)
{
	bool is_const_int =((strcmp(type,"sp") == 0) || (strcmp(type,"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(type,"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"lOpen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"l1pen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"wav") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"blkdiag") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin_test") == 0)));
	
	return is_const_int;
}


bool isConstraintNameReal(const char * type)
{
	return ((strcmp(type,"normcol") == 0) || (strcmp(type,"normlin")==0));
}

bool isConstraintNameMat(const char * type)
{
	return ((strcmp(type,"supp") == 0) || (strcmp(type,"const")==0));
}

faust_constraint_name getEquivalentConstraint(const char * type)
{	
	
	if (!strcmp(type,"sp"))	
		return CONSTRAINT_NAME_SP;
	if (!strcmp(type,"spcol"))	
		return CONSTRAINT_NAME_SPCOL;
	if  (!strcmp(type,"splin"))	
		return CONSTRAINT_NAME_SPLIN;
	if  (!strcmp(type,"normcol"))	
		return CONSTRAINT_NAME_NORMCOL;
	if  (!strcmp(type,"splincol"))	
		return CONSTRAINT_NAME_SPLINCOL;
	if  (!strcmp(type,"l0pen"))	
		return CONSTRAINT_NAME_L0PEN;
	if  (!strcmp(type,"l1pen"))	
		return CONSTRAINT_NAME_L1PEN;
	if  (!strcmp(type,"const"))	
		return CONSTRAINT_NAME_CONST;
	if  (!strcmp(type,"wav"))	
		return CONSTRAINT_NAME_WAV;
	if  (!strcmp(type,"sppos"))	
		return CONSTRAINT_NAME_SP_POS;
	if  (!strcmp(type,"blkdiag"))	
		return CONSTRAINT_NAME_BLKDIAG;
	if  (!strcmp(type,"splin_test"))	
		return CONSTRAINT_NAME_SPLIN_TEST;
	if  (!strcmp(type,"supp"))	
		return CONSTRAINT_NAME_SUPP;
	if  (!strcmp(type,"normlin"))	
		return CONSTRAINT_NAME_NORMLIN;
	if  (!strcmp(type,"toeplitz"))	
		return CONSTRAINT_NAME_TOEPLITZ;

	
	handleError("faust_constraint_generic","getEquivalentConstraint : Unknown type of constraint");
}


int getTypeConstraint(const char * type)
{
	bool is_const_int = isConstraintNameInt(type);
	bool is_const_real = isConstraintNameReal(type);
	bool is_const_mat = isConstraintNameMat(type);
				

	if (is_const_int)
		return 0;
	if (is_const_real)
		return 1;
	if (is_const_mat)
		return 2;
	handleError("faust_constraint_generic","::add_constraint : invalid constraint type");
}

