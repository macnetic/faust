#ifndef __FAUST_INIT_FROM_MATIO_PARAMS_H__
#define __FAUST_INIT_FROM_MATIO_PARAMS_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

class faust_params_palm;
class faust_params;
class faust_constraint_generic;

void init_params_palm_from_matiofile(faust_params_palm& params, const char* fileName, const char* variableName);
void init_params_from_matiofile(faust_params& params, const char* fileName, const char* variableName);

void add_constraint(std::vector<const faust_constraint_generic*> & consS,matvar_t* cons_var);
void Display_params(faust_params & params);
#endif
