#ifndef __FAUST_INIT_FROM_MATIO_PARAMS_H__
#define __FAUST_INIT_FROM_MATIO_PARAMS_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

#include "faust_Params.h"
#include "faust_ParamsPalm.h"



template<typename FPP,Device DEVICE> class ParamsPalm;
template<typename FPP,Device DEVICE> class Params;
template<typename FPP,Device DEVICE> class ConstraintGeneric;

template<typename FPP,Device DEVICE>
void init_params_palm_from_matiofile(Faust::ParamsPalm<FPP,DEVICE>& params, const char* fileName, const char* variableName);

/** \brief load data matrix from ".mat file"
 * \param params
 * \param fileName
 * \param variableName
 */
template<typename FPP,Device DEVICE>
void init_params_from_matiofile(Faust::Params<FPP,DEVICE>& params, const char* fileName, const char* variableName);

template<typename FPP,Device DEVICE>
void add_constraint(std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> & consS,matvar_t* cons_var);

template<typename FPP,Device DEVICE>
void Display_params(Faust::Params<FPP,DEVICE> & params);

#include "faust_init_from_matio_params.hpp"
#endif
