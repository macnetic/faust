#ifndef LINALGEBRA_H
#define LINALGEBRA_H

#include "faust_constant.h"
#include "faust_spmat.h"

class faust_mat;
class faust_vec;
//class faust_spmat;
class faust_core;

#endif


//////////FONCTION faust_mat - faust_mat ////////////////////

// C = A * B; 
//l'objet C doit etre different de A et B
 void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C);
 
 // C = A + B
 void add(const faust_mat & A, const faust_mat & B, faust_mat & C);
 
 // C = alpha * A*B + beta * C;
 // l'objet C doit etre different de A et B
 //void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta);
 
 
faust_vec solve(const faust_mat & A, const faust_vec & v);
void solve(const faust_spmat & A,faust_vec & x, const faust_vec & y);
 
 
 // C = alpha *op(A)*op(B) + beta * C;
 // op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'
 // op(B) = B si typeB='N', op(B) = transpose(B) si typeB='T'
 // l'objet C doit etre different de A et B
 void gemm(const faust_mat & A,const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA, char  typeB);
 void gemv(const faust_mat & A,const faust_vec & x,faust_vec & y,const faust_real & alpha, const faust_real & beta, char typeA);
 
 faust_real power_iteration(const faust_mat & A, const int nbr_iter_max,faust_real threshold,int & flag);
 faust_real power_iteration(const  faust_core & A, const int nbr_iter_max,faust_real threshold, int & flag);
 // C = op(A) * B if typeMult = 'R'
 // C = B * op(A) sinon
 // op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'
 void multiply(const faust_core & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char typeA, char typeMult);

// non-member operators declarations
faust_vec operator*(const faust_core& f, const faust_vec& v);
faust_mat operator*(const faust_core& f, const faust_mat& M);
faust_vec operator*(const faust_mat& M, const faust_vec& v);
