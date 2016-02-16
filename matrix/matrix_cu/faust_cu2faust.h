#ifndef __FAUST_CU2FAUST_H__
#define __FAUST_CU2FAUST_H__
#include "faust_cuda.h"
template<typename faust_real> class faust_vec;
template<typename faust_real> class faust_mat;
template<typename faust_real> class faust_cu_vec;
#include "faust_cu_mat.h"

#ifdef __COMPILE_SPMAT__
template<typename faust_real> class faust_spmat;
template<typename faust_real> class faust_cu_spmat;
#endif

template<typename T, typename U>
void faust_cu2faust(faust_vec<T>& v, const faust_cu_vec<U>& cu_v, cudaStream_t stream=0);
template<typename T, typename U>
void faust_cu2faust(faust_mat<T>& M, const faust_cu_mat<U>& cu_M, cudaStream_t stream=0);
#ifdef __COMPILE_SPMAT__
template<typename T, typename U>
void faust_cu2faust(faust_spmat<T>& S, const faust_cu_spmat<U>& cu_S, cudaStream_t stream=0);
#endif


#include "faust_cu2faust.hpp"

#endif
