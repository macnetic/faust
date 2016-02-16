#ifndef __FAUST_CORE_ALGEBRA_CU_H__
#define __FAUST_CORE_ALGEBRA_CU_H__

#include "faust_constant.h"
#include "faust_spmat.h"

template<typename T> class faust_cu_mat;
template<typename T> class faust_cu_vec;
template<typename T> class faust_cu_core;


template<typename T> 
T power_iteration(const faust_core_cu<T> & A, const int nbr_iter_max, const T threshold, int & flag);




template<typename T> 
faust_cu_vec<T> operator*(const faust_core_cu<T>& f, const faust_cu_vec<T> & v);

template<typename T> 
faust_cu_mat<T> operator*(const faust_core_cu<T>& f, const faust_cu_mat<T> & M);

#include "faust_core_algebra_cu.hpp"

#endif
