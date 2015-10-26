#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"


class faust_mat_generic
{
public:	

faust_unsigned_int getNbRow() const {return dim1;}
faust_unsigned_int getNbCol() const {return dim2;}

faust_mat_generic() : dim1(0), dim2(0) {}
faust_mat_generic(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_){} 

protected:
	faust_unsigned_int dim1;
    faust_unsigned_int dim2;
	
	
};

#endif