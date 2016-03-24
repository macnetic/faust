#ifndef __FAUST_LINEAR_OPERATOR_H__
#define __FAUST_LINEAR_OPERATOR_H__

#include "faust_constant.h"
#include "faust_mat.h"


template<typename T> class faust_mat;

template<typename FPP>
class faust_linear_operator
{
	public:	
	
	virtual faust_unsigned_int getNbRow() const=0;
	virtual faust_unsigned_int getNbCol() const=0; 
	virtual void transpose()=0;
	virtual void mult(faust_mat<FPP> Y,const  faust_mat<FPP> X)const=0;
	


	
};
#endif
