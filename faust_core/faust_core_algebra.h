#ifndef FAUST_CORE_ALGEBRA_H
#define FAUST_CORE_ALGEBRA_H

#include "faust_constant.h"
#include "faust_spmat.h"

class faust_mat;
class faust_vec;
//class faust_spmat;
class faust_core;

#endif

faust_real power_iteration(const  faust_core & A, const int nbr_iter_max,faust_real threshold, int & flag);


 void multiply(const faust_core & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char typeA, char typeMult);
 void faust_solve(const faust_core & A,faust_vec & x, const faust_vec & y);



faust_vec operator*(const faust_core& f, const faust_vec& v);
faust_mat operator*(const faust_core& f, const faust_mat& M);