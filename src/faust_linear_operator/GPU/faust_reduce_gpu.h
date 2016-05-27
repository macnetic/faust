
#ifndef __FAUST_CU_REDUCE_H__
#define __FAUST_CU_REDUCE_H__

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>

/*template<typename faust_real> class faust_cu_vec;
template<typename faust_real> class faust_cu_matrix;

template<typename faust_real>
faust_real faust_cu_reduce(const faust_cu_vec<faust_real>& v);

template<typename faust_real>
faust_real faust_cu_reduce(const faust_cu_matrix<faust_real>& M);*/


template<typename FPP>
FPP faust_cu_sum(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_max(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_min(const FPP* data, const int nb_el);

template<typename FPP>
FPP faust_cu_norm(const FPP* data, const int nb_el);



#endif
