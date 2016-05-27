#ifndef __FAUST_INIT_FROM_MATIO_CORE_H__
#define __FAUST_INIT_FROM_MATIO_CORE_H__

#include "faust_init_from_matio.h"
#include "matio.h"
#include <vector>

// class Faust::Transform;
// class Faust::MatDense;
template<typename FPP,Device DEVICE> class Transform;
template<typename FPP,Device DEVICE> class MatDense;


template<typename FPP,Device DEVICE>
void init_faust_core_from_matiofile(Faust::Transform<FPP,DEVICE>& core, const char* fileName, const char* variableName);
template<typename FPP,Device DEVICE>
void init_faust_core_from_matvar(Faust::Transform<FPP,DEVICE>& core, matvar_t* cell_var );
template<typename FPP,Device DEVICE>
void init_faust_data_from_matiofile(std::vector<Faust::MatDense<FPP,DEVICE> >& full_mat, std::vector<Faust::Transform<FPP,DEVICE> >& core, const char* fileName, const char* variableName);

#ifdef COMPILE_GPU
template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Gpu> core, const char* fileName, const char* variableName);
#endif


//template<typename FPP,Device DEVICE>
//void write_faust_spmat_list_into_matfile(const faust_spmat<FPP,DEVICE>& M, const char* fileName, const char* variableName);

template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Cpu> core, const char* fileName, const char* variableName);

template<typename FPP>
void write_faust_core_into_matvar(const Faust::Transform<FPP,Cpu> core, matvar_t** matvar, const char* variableName);





#include "faust_init_from_matio_core.hpp"

#endif
