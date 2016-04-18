#ifndef __PROX_CU_H__
#define __PROX_CU_H__

#include "faust_cu_mat.h"
#include "faust_constant.h"
#include "faust_exception.h"
#include "faust_vec.h"

template<typename T>
void prox_sp(faust_cu_mat<T> & M,faust_unsigned_int k);
template<typename T>
void prox_sp_pos(faust_cu_mat<T> & M,faust_unsigned_int k);
template<typename T>
void prox_spcol(faust_cu_mat<T> & M,faust_unsigned_int k);
template<typename T>
void prox_splin(faust_cu_mat<T> & M,faust_unsigned_int k);
template<typename T>
void prox_splincol(faust_cu_mat<T> & M,faust_unsigned_int k);


template<typename T>
void prox_normcol(faust_cu_mat<T> & M,T s);
template<typename T>
void prox_normlin(faust_cu_mat<T> & M,T s);
template<typename T>
void prox_supp(faust_cu_mat<T> & M, const faust_cu_mat<T> & supp);
template<typename T>
void prox_blkdiag(faust_cu_mat<T> & M,faust_unsigned_int k);
template<typename T>
void prox_splincol(faust_cu_mat<T> & M,faust_unsigned_int k);




#include "prox_cu.hpp"

#endif
