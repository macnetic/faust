#ifndef __FAUST_INIT_FROM_MATIO_H__
#define __FAUST_INIT_FROM_MATIO_H__

#include "matio.h"
#include <vector>
#include "faust_params_palm.h"
#include "faust_params.h"


class faust_mat;

matvar_t* faust_matio_read_variable(mat_t* file, const char* variableName);
void init_faust_mat_from_matio_mat(faust_mat& M, const char* fileName, const char* variableName);
double init_faust_mat_from_matio_double(const char* fileName, const char* variableName);
int init_faust_mat_from_matio_int(const char* fileName, const char* variableName);
bool init_faust_mat_from_matio_bool(const char* fileName, const char* variableName);
void init_faust_mat_vector_from_matiofile( std::vector<faust_mat> & vec_M, const char* fileName, const char* variableName);
void init_params_palm_from_matiofile(faust_params_palm& params, const char* fileName, const char* variableName);
void init_params_from_matiofile(faust_params& params, const char* fileName, const char* variableName);

void write_faust_mat_into_matfile( faust_mat& M, const char* fileName, const char* variableName);

#endif
