#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"
#include "faust_linear_operator.h"
/**
 * \class faust_mat_generic faust_mat_generic.h
 * \brief This faust_mat_generic class serves as a base class for the derived class faust_mat and faust_spmat .
 * 	Member variable dim1 and dim2 correspond to the dimension of the matrix
 *
*/

template<typename FPP>
class faust_mat_generic : public faust_linear_operator<FPP>
{
	public:	
	faust_mat_generic() : dim1(0), dim2(0) {}
	faust_mat_generic(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_){}

	faust_unsigned_int getNbRow() const {return dim1;}
	faust_unsigned_int getNbCol() const {return dim2;} 
	
	void resize(const faust_unsigned_int dim1_,const faust_unsigned_int dim2_){dim1=dim1_;dim2=dim2_;}

	protected:
	faust_unsigned_int dim1;
   	faust_unsigned_int dim2;	
};
#endif
