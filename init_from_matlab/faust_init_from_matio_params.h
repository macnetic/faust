#ifndef __FAUST_INIT_FROM_MATIO_PARAMS_H__
#define __FAUST_INIT_FROM_MATIO_PARAMS_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>
#include "faust_params.h"
#include "faust_params_palm.h"


template<typename T> class faust_params_palm;
template<typename T> class faust_params;
class faust_constraint_generic;

template<typename T>
void init_params_palm_from_matiofile(faust_params_palm<T>& params, const char* fileName, const char* variableName);

/** \brief load data matrix from ".mat file"
 * \param params
 * \param fileName
 * \param variableName
 */
template<typename T>
void init_params_from_matiofile(faust_params<T>& params, const char* fileName, const char* variableName);

template<typename T>
void add_constraint(std::vector<const faust_constraint_generic*> & consS,matvar_t* cons_var);

template<typename T>
void Display_params(faust_params<T> & params);

#include "faust_init_from_matio_params.hpp"
#endif
