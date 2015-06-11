#ifndef LINALGEBRA_H
#define LINALGEBRA_H
#include "faust_mat.h"


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
 
 
 // C = alpha *op(A)*op(B) + beta * C;
 // op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'
 // op(B) = B si typeB='N', op(B) = transpose(B) si typeB='T'
 // l'objet C doit etre different de A et B
 void gemm(faust_mat & A, faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA, char  typeB);


