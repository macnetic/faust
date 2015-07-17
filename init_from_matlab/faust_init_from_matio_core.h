#ifndef __FAUST_INIT_FROM_MATIO_CORE_H__
#define __FAUST_INIT_FROM_MATIO_CORE_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

class faust_core;
class faust_mat;

void init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName);
void init_faust_core_from_matvar(faust_core& core, matvar_t* cell_var );
void init_faust_data_from_matiofile(std::vector<faust_mat>& full_mat, std::vector<faust_core>& core, const char* fileName, const char* variableName);


#endif
