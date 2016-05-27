#ifndef __FAUST_Transform_ALGEBRA_GPU_H__
#define __FAUST_Transform_ALGEBRA_GPU_H__

#include "faust_constant.h"

template<typename FPP,Device DEVICE> class MatDense;
template<typename FPP,Device DEVICE> class MatSparse;
template<typename FPP,Device DEVICE> class Vect;
template<typename FPP,Device DEVICE> class Transform;

// power iteration was not tasted with Faust::Transform<FPP,Gpu> because it seems that it is not tested
//template<typename FPP>
//FPP power_iteration(const Faust::Transform<FPP,Gpu> & A, const int nbr_iter_max, const FPP threshold, int & flag, cublasHandle_t cublasHandle);



//currently operator * (multiplication) can not be overloaded because GPU multiplication implies taking extra arguments cublashandle (dense matrix)  and cusparsehandle (sparse matrix)
//whereas operator must take either one or two arguments

//template<typename T>
//faust_cu_vec<T> operator*(const Faust::Transform_cu<T>& f, const faust_cu_vec<T> & v);

//template<typename T>
//faust_cu_mat<T> operator*(const Faust::Transform_cu<T>& f, const faust_cu_mat<T> & M);

#include "faust_transform_algebra_gpu.hpp"

#endif
