
#include <vector>
#include <iostream>
#include "algorithm"
#include "faust_cu2faust.h"
#include "faust_mat.h"

// const char * interface_prox_name="prox : ";

template<typename T>
inline bool partial_sort_comp (const std::pair<int, T>& pair1, const std::pair<int, T>& pair2) 
{ 
   return fabs(pair1.second) > fabs(pair2.second); 
}

template<typename T>
void sort_idx(const std::vector<T> &v, std::vector<int>& idx, int s) 
{
      std::vector<std::pair<int, T> > vec_pair(v.size());
      for (int i=0 ; i<v.size() ; i++)
          vec_pair[i] = std::make_pair(i,v[i]);
     
      std::partial_sort(vec_pair.begin(), vec_pair.begin()+s, vec_pair.end(),partial_sort_comp<T>);
      idx.resize(s);
      for (int i=0 ; i<s ; i++)
          idx[i]=vec_pair[i].first;
}

template<typename T>
void prox_sp(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_sp(M,k)
   cu_M = M;
}

template<typename T>
void prox_spcol(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_spcol(M,k)
   cu_M = M;
}

template<typename T>
void prox_splin(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_splin(M,k)
   cu_M = M;
}

template<typename T>
void prox_normcol(faust_cu_mat<T>& cu_M, T s)
{
	faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_normcol(M,s)
   cu_M = M;
}

template<typename T>
void prox_normlin(faust_cu_mat<T> & cu_M, T s)
{
	faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_normlin(M,s)
   cu_M = M;	
}

template<typename T>
void prox_sp_pos(faust_cu_mat<T>& cu_M, faust_unsigned_int k)
{
   faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_sp_pos(M,k)
   cu_M = M;
}


// M needs to be square and k must divide dimension of M
template<typename T>
void prox_blkdiag(faust_cu_mat<T>& cu_M, int k)
{   
   faust_mat M;
   faust_cu2faust(M, cu_M);
   prox_blkdiag(M,k)
   cu_M = M;	
}

template<typename T>
void prox_supp(faust_cu_mat<T>& cu_M, const faust_cu_mat<T>& cu_supp)
{   
   faust_mat M,supp;
   faust_cu2faust(M, cu_M);
   faust_cu2faust(supp, cu_supp);
   prox_supp(M,supp)
   cu_M = M;	
}

