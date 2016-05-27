#ifndef __PROX_CU_H__
#define __PROX_CU_H__

#include "faust_MatDense_gpu.h"
#include "faust_constant.h"
#include "faust_exception.h"


template<typename FPP>
void prox_sp(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_sp_pos(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_spcol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splin(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splincol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);


template<typename FPP>
void prox_normcol(Faust::MatDense<FPP,Gpu> & M,FPP s);
template<typename FPP>
void prox_normlin(Faust::MatDense<FPP,Gpu> & M,FPP s);
template<typename FPP>
void prox_supp(Faust::MatDense<FPP,Gpu> & M, const Faust::MatDense<FPP,Gpu> & supp);
template<typename FPP>
void prox_blkdiag(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splincol(Faust::MatDense<FPP,Gpu> & M,faust_unsigned_int k);




#include "prox_gpu.hpp"

#endif
