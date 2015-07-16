#ifndef __FAUST_INIT_FROM_MATIO_CORE_H__
#define __FAUST_INIT_FROM_MATIO_CORE_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

class faust_core;

void init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName);

#endif
