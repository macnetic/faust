#ifndef SPBLASHANDLE_CPU_H
#define SPBLASHANDLE_CPU_H
#include "faust_constant.h"


template <Device DEVICE> class SpBlasHandle;
 

//BlasHandle<Cpu> is useless,
// its unique purpose is for code factorisation with GPU code

template <>
class SpBlasHandle<Cpu>
{

public :

	SpBlasHandle(){} 
 	SpBlasHandle(const SpBlasHandle<Cpu> & blashandle){}
	void operator=(SpBlasHandle<Cpu> const& blashandle){};
	void Init(){}
	void Destroy(){}

	
	
};



#endif
