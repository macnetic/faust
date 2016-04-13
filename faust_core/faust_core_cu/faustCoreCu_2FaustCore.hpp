#ifndef __FAUSTCORECU_2FAUSTCORE_HPP__
#define __FAUSTCORECU_2FAUSTCORE_HPP__

#include "faust_core.h"
#include "faust_core_cu.h"
#include "faust_cu2faust.h"


template<typename T, typename U>
void faust_cu2faust(faust_core<T>& fcore, const faust_core_cu<U>& cu_fcore)
{
	
	std::vector< faust_spmat<T> > list_spmat(cu_fcore.size());
	std::vector< faust_cu_spmat<U> > list_cu_spmat;
	cu_fcore.get_facts(list_cu_spmat);

	for (int i=0;i<cu_fcore.size();i++)
	{
		faust_cu2faust(list_spmat[i],list_cu_spmat[i]);
	}
	faust_core<T> fcorebis(list_spmat);
	fcore=fcorebis;	
}

#endif

