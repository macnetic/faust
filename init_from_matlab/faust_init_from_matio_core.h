#ifndef __FAUST_INIT_FROM_MATIO_CORE_H__
#define __FAUST_INIT_FROM_MATIO_CORE_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

// class faust_core;
// class faust_mat;
template<typename T> class faust_core;
template<typename T> class faust_mat;

template<typename T>
void init_faust_core_from_matiofile(faust_core<T>& core, const char* fileName, const char* variableName);
template<typename T>
void init_faust_core_from_matvar(faust_core<T>& core, matvar_t* cell_var );
template<typename T>
void init_faust_data_from_matiofile(std::vector<faust_mat<T> >& full_mat, std::vector<faust_core<T> >& core, const char* fileName, const char* variableName);

#include "faust_init_from_matio_core.hpp"

#endif
