#ifndef PROX_H
#define PROX_H
#include "faust_mat.h"


void prox_sp(faust_mat & M,int k);
void prox_sp_pos(faust_mat & M,int k);
void prox_spcol(faust_mat & M,int k);
void prox_splin(faust_mat & M,int k);
void prox_normcol(faust_mat & M, faust_real s);

#endif
