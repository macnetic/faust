#ifndef __FAUST_INIT_FROM_MATIO_MAT_H__
#define __FAUST_INIT_FROM_MATIO_MAT_H__

#include "matio.h"
#include <vector>


class faust_mat;

void init_faust_mat_from_matio(faust_mat& M, const char* fileName, const char* variableName);
void init_mat_from_matvar(faust_mat & M,matvar_t* var);
void write_faust_mat_into_matfile( faust_mat& M, const char* fileName, const char* variableName);
void init_faust_mat_vector_from_matiofile( std::vector<faust_mat> & vec_M, const char* fileName, const char* variableName);
#endif
