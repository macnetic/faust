#ifndef PROX_H
#define PROX_H

#include "faust_MatDense.h"
#include "faust_constant.h"
#include "faust_exception.h"
#include "faust_Vect.h"




template<typename FPP>
bool partial_sort_comp (const std::pair<int, FPP>& pair1, const std::pair<int, FPP>& pair2);

template<typename FPP>
void sort_idx(const std::vector<FPP> &v, std::vector<int>& idx, int s);

template<typename FPP>
void prox_sp(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_sp_pos(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_spcol(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splin(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_spcol_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splin_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);
template<typename FPP>
void prox_splincol(Faust::MatDense<FPP,Cpu> &M,faust_unsigned_int k);
template<typename FPP>
void prox_normcol(Faust::MatDense<FPP,Cpu> & M,FPP s);
template<typename FPP>
void prox_normlin(Faust::MatDense<FPP,Cpu> & M,FPP s);
template<typename FPP>
void prox_supp(Faust::MatDense<FPP,Cpu> & M, const Faust::MatDense<FPP,Cpu> & supp);
template<typename FPP>
void prox_blkdiag(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k);



#include "prox.hpp"

#endif
