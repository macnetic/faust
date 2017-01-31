#ifndef __FAUST_LINEAR_OPERATOR_HPP__
#define __FAUST_LINEAR_OPERATOR_HPP__

#include <complex>
#include <typeinfo>
#include "faust_exception.h"


template<typename FPP,Device DEVICE>
bool Faust::LinearOperator<FPP,DEVICE>::isReal() const
{
	bool isReal = (typeid(FPP) == typeid(double));
	     	isReal = (isReal || (typeid(FPP) == typeid(float)) );

		bool isComplex = (typeid(FPP) == typeid(std::complex<double>));
	     	isComplex = (isComplex || (typeid(FPP) == typeid(std::complex<float>)) );	
	
		if  ( (!isComplex) && (!isReal) )
		{
			handleError("linearOperator","isReal : unknown type of scalar");
		}
	
		return isReal;
	
}

#endif
