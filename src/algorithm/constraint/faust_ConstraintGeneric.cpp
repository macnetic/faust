#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include <typeinfo>
#include "faust_exception.h"

#include "faust_ConstraintMat.h"

// modif AL AL
#include "faust_ConstraintType.h"




bool is_constraint_name_int(const char * type)
{
	bool is_const_int =((strcmp(type,"sp") == 0) || (strcmp(type,"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(type,"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"blkdiag") == 0)));


	return is_const_int;
}


bool is_constraint_name_real(const char * type)
{
	return ((strcmp(type,"normcol") == 0) || (strcmp(type,"normlin")==0));
}

bool is_constraint_name_mat(const char * type)
{
	return ((strcmp(type,"supp") == 0) || (strcmp(type,"const")==0));
}

faust_constraint_name get_equivalent_constraint(const char * type)
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
	if  (!strcmp(type,"const"))
		return CONSTRAINT_NAME_CONST;
	if  (!strcmp(type,"sppos"))
		return CONSTRAINT_NAME_SP_POS;
	if  (!strcmp(type,"blkdiag"))
		return CONSTRAINT_NAME_BLKDIAG;
	if  (!strcmp(type,"supp"))
		return CONSTRAINT_NAME_SUPP;
	if  (!strcmp(type,"normlin"))
		return CONSTRAINT_NAME_NORMLIN;


	handleError("Faust::ConstraintGeneric","get_equivalent_constraint : Unknown type of constraint");
}


int get_type_constraint(const char * type)
{
	bool is_const_int = is_constraint_name_int(type);
	bool is_const_real = is_constraint_name_real(type);
	bool is_const_mat = is_constraint_name_mat(type);


	if (is_const_int)
		return 0;
	if (is_const_real)
		return 1;
	if (is_const_mat)
		return 2;
	handleError("Faust::ConstraintGeneric","::add_constraint : invalid constraint type");
}

