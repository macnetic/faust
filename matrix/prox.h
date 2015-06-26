#ifndef PROX_H
#define PROX_H
#include "faust_mat.h"

//void new_prox_sp(faust_mat & M,int k);


void prox_sp(faust_mat & M,int k);
void prox_sp_pos(faust_mat & M,int k);
void prox_spcol(faust_mat & M,int k);
void prox_splin(faust_mat & M,int k);
void prox_normcol(faust_mat & M,faust_real s);
void prox_normlin(faust_mat & M,faust_real s);
void prox_supp(faust_mat & M, const faust_mat & supp);
void prox_blkdiag(faust_mat & M,int k);
void prox_toeplitz(faust_mat & M, int k);


void prox_spcol_old(faust_mat & M,int k);
void prox_splin_old(faust_mat & M,int k);
void prox_sp_old(faust_mat & M,int k);
void prox_sp_old_old(faust_mat & M,int k);

#endif
