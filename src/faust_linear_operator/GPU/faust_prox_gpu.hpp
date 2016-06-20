#ifndef __FAUST_PROX_CU_HPP__
#define __FAUST_PROX_CU_HPP__

#include <vector>
#include <iostream>
#include "algorithm"
#include "faust_gpu2cpu.h"
#include "faust_MatDense.h"
#include "faust_prox.h"

// const char * interface_prox_name="prox : ";


template<typename FPP>
void Faust::prox_sp(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_sp(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_spcol(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_spcol(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_splin(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_splin(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_normcol(Faust::MatDense<FPP,Gpu>& cu_M, FPP s)
{
    Faust::MatDense<FPP,Cpu> M;
    faust_gpu2cpu(M, cu_M);
    Faust::prox_normcol(M,s);
    cu_M = M;
}

template<typename FPP>
void Faust::prox_normlin(Faust::MatDense<FPP,Gpu> & cu_M, FPP s)
{
	Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_normlin(M,s);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_sp_pos(Faust::MatDense<FPP,Gpu>& cu_M, faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_sp_pos(M,k);
   cu_M = M;
}


// M needs to be square and k must divide dimension of M
template<typename FPP>
void Faust::prox_blkdiag(Faust::MatDense<FPP,Gpu>& cu_M, int k)
{
   Faust::MatDense<FPP,Cpu> M;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_blkdiag(M,k);
   cu_M = M;
}

template<typename FPP>
void Faust::prox_supp(Faust::MatDense<FPP,Gpu>& cu_M, const Faust::MatDense<FPP,Gpu>& cu_supp)
{
   Faust::MatDense<FPP,Cpu> M,supp;
   faust_gpu2cpu(M, cu_M);
   faust_gpu2cpu(supp, cu_supp);
   Faust::prox_supp(M,supp);
   cu_M = M;
}



template<typename FPP>
void Faust::prox_splincol(Faust::MatDense<FPP,Gpu>& cu_M, const faust_unsigned_int k)
{
   Faust::MatDense<FPP,Cpu> M,supp;
   faust_gpu2cpu(M, cu_M);
   Faust::prox_splincol(M,k);
   cu_M = M;
}
#endif


