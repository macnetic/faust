#ifndef FAUST_CORE_ALGEBRA_H
#define FAUST_CORE_ALGEBRA_H

#include "faust_constant.h"
#include "faust_spmat.h"

template<typename T> class faust_mat;
template<typename T> class faust_vec;
template<typename T> class faust_core;


template<typename T> 
T power_iteration(const faust_core<T> & A, const int nbr_iter_max,T threshold, int & flag);




template<typename T> 
faust_vec<T> operator*(const faust_core<T>& f, const faust_vec<T> & v);

template<typename T> 
faust_mat<T> operator*(const faust_core<T>& f, const faust_mat<T> & M);

#include "faust_core_algebra.hpp"

#endif
