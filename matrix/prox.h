#ifndef PROX_H
#define PROX_H
#include "faust_mat.h"
#include "faust_constant.h"

//void new_prox_sp(faust_mat & M,int k);


bool partial_sort_comp (const std::pair<int, faust_real>& pair1, const std::pair<int, faust_real>& pair2);

void sort_idx(const std::vector<faust_real> &v, std::vector<int>& idx, int s); 

void prox_sp(faust_mat & M,faust_unsigned_int k);
void prox_sp_pos(faust_mat & M,faust_unsigned_int k);
void prox_spcol(faust_mat & M,faust_unsigned_int k);
void prox_splin(faust_mat & M,faust_unsigned_int k);
void prox_normcol(faust_mat & M,faust_real s);
void prox_normlin(faust_mat & M,faust_real s);
void prox_supp(faust_mat & M, const faust_mat & supp);
void prox_blkdiag(faust_mat & M,faust_unsigned_int k);
//void prox_toeplitz(faust_mat & M, int k);

void prox_sp_normfree(faust_mat & M,faust_unsigned_int k);
void prox_sp_pos_normfree(faust_mat & M,faust_unsigned_int k);
void prox_spcol_normfree(faust_mat & M,faust_unsigned_int k);
void prox_splin_normfree(faust_mat & M,faust_unsigned_int k);
void prox_supp_normfree(faust_mat & M,const faust_mat & supp);



#endif
