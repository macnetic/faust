#ifndef __FAUST_INIT_FROM_MATIO_MAT_H__
#define __FAUST_INIT_FROM_MATIO_MAT_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include "faust_constant.h"
#include <vector>
#include "faust_mat.h"
#include "faust_spmat.h" 

template<typename T> class faust_spmat;
template<typename T> class faust_mat;


template<typename T>
void init_faust_mat_from_matio(faust_mat<T>& M, const char* fileName, const char* variableName);
template<typename T>
void init_faust_spmat_from_matio(faust_spmat<T>& S, const char* fileName, const char* variableName);
template<typename T>
void write_faust_mat_list_into_matfile(const std::vector< faust_mat<T> >& M, const char* fileName, const char* variableName);
template<typename T>
void init_faust_mat_vector_from_matiofile( std::vector<faust_mat<T> > & vec_M, const char* fileName, const char* variableName);
template<typename T>
void init_mat_from_matvar(faust_mat<T> & M,matvar_t** var);
template<typename T>
void init_spmat_from_matvar(faust_spmat<T> & M,matvar_t* var);






template<typename T>
void write_faust_mat_list_into_matvar(const std::vector<faust_mat<T> >& M,matvar_t** matvar, const char* variableName);


//passer l adresse du pointeur en parametre, pas un pointeur de pointeur pour matvar_t** matvar
template<typename T>
void write_faust_mat_into_matvar(const faust_mat<T>& M,matvar_t** matvar, const char* variableName);
template<typename T>
void write_faust_mat_into_matfile(const faust_mat<T>& M, const char* fileName, const char* variableName);

/// A MODIFIER : pour l instant les matrices creuse sont stockés en matrices dense dans matlab
template<typename T>
void write_faust_spmat_list_into_matfile(const faust_spmat<T>& M, const char* fileName, const char* variableName);
template<typename T>
void write_faust_spmat_list_into_matvar(const std::vector<faust_spmat<T> >& M,matvar_t** matvar, const char* variableName);





#include "faust_init_from_matio_mat.hpp"

#endif
