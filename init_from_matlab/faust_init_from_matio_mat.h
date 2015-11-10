#ifndef __FAUST_INIT_FROM_MATIO_MAT_H__
#define __FAUST_INIT_FROM_MATIO_MAT_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include "faust_constant.h"
#include <vector>

template<typename T> class faust_spmat;
template<typename T> class faust_mat;

template<typename T>
void init_faust_mat_from_matio(faust_mat<T>& M, const char* fileName, const char* variableName);
template<typename T>
void init_faust_spmat_from_matio(faust_spmat<T>& S, const char* fileName, const char* variableName);
template<typename T>
void write_faust_mat_into_matfile( faust_mat& M<T>, const char* fileName, const char* variableName);
template<typename T>
void init_faust_mat_vector_from_matiofile( std::vector<faust_mat<T> > & vec_M, const char* fileName, const char* variableName);
template<typename T>
void init_mat_from_matvar(faust_mat & M<T>,matvar_t* var);
template<typename T>
void init_spmat_from_matvar(faust_spmat<T> & M,matvar_t* var);


#include "faust_init_from_matio_mat.hpp"

#endif
