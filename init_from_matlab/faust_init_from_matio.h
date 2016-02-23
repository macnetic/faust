#ifndef __FAUST_INIT_FROM_MATIO_H__
#define __FAUST_INIT_FROM_MATIO_H__

#include "matio.h"


matvar_t* faust_matio_read_variable(const char* fileName, const char* variableName);
double init_double_from_matio(const char* fileName, const char* variableName);
int init_int_from_matio(const char* fileName, const char* variableName);
bool init_bool_from_matio(const char* fileName, const char* variableName);


#endif
