#ifndef BLASHANDLE_CPU_H
#define BLASHANDLE_CPU_H
#include "faust_constant.h"


template <Device DEVICE> class BlasHandle;
 

//BlasHandle<Cpu> is useless,
// its unique purpose is for code factorisation with GPU code

template <>
class BlasHandle<Cpu>
{

public :

	BlasHandle(){} 
 	BlasHandle(const BlasHandle<Cpu> & blashandle){}
	void operator=(BlasHandle<Cpu> const& blashandle){};
	void Init(){}
	void Destroy(){}

	
	
};



#endif
