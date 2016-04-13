#ifndef __FAUSTCORECU_2FAUSTCOREH__
#define __FAUSTCORECU_2FAUSTCOREH__


template<typename T> class faust_core_cu;
template<typename T> class faust_core;

template<typename T, typename U>
void faust_cu2faust(faust_core<T>& fcore, const faust_core_cu<U>& cu_fcore);

#include "faustCoreCu_2FaustCore.hpp"
#endif
