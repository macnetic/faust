#ifndef __FAUST_PROX_CU_HPP__
#define __FAUST_PROX_CU_HPP__

#include <vector>
#include <iostream>
#include "algorithm"
#include "faust_cu2faust.h"
#include "faust_mat.h"
#include "prox.h"

// const char * interface_prox_name="prox : ";


template<typename T>
void prox_sp(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_sp(M,k);
   cu_M = M;
}

template<typename T>
void prox_spcol(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_spcol(M,k);
   cu_M = M;
}

template<typename T>
void prox_splin(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_splin(M,k);
   cu_M = M;
}

template<typename T>
void prox_normcol(faust_cu_mat<T>& cu_M, T s)
{
	faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_normcol(M,s);
   cu_M = M;
}

template<typename T>
void prox_normlin(faust_cu_mat<T> & cu_M, T s)
{
	faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_normlin(M,s);
   cu_M = M;	
}

template<typename T>
void prox_sp_pos(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_sp_pos(M,k);
   cu_M = M;
}


// M needs to be square and k must divide dimension of M
template<typename T>
void prox_blkdiag(faust_cu_mat<T>& cu_M, int k)
{   
   faust_mat<T> M;
   faust_cu2faust(M, cu_M);
   prox_blkdiag(M,k);
   cu_M = M;	
}

template<typename T>
void prox_supp(faust_cu_mat<T>& cu_M, const faust_cu_mat<T>& cu_supp)
{   
   faust_mat<T> M,supp;
   faust_cu2faust(M, cu_M);
   faust_cu2faust(supp, cu_supp);
   prox_supp(M,supp);
   cu_M = M;	
}

#endif
