#ifndef LINALGEBRA_H
#define LINALGEBRA_H

#include "faust_constant.h"
#include "faust_spmat.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

template<typename T> class faust_vec;
template<typename T> class faust_mat;
template<typename T> class faust_spmat;





//////////FONCTION faust_mat<T> - faust_mat<T> ////////////////////

// C = A * B; 
//l'objet C doit etre different de A et B
template<typename T>
 void multiply(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
 
 // C = A + B
// template<typename T> 
 // void add(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
 
 // C = alpha * A*B + beta * C;
 // l'objet C doit etre different de A et B
 //void gemm(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta);
 
template<typename T> 
faust_vec<T> solve(const faust_mat<T> & A, const faust_vec<T> & v);

// template<typename T> 
// void sp_solve(const faust_spmat<T> & A,faust_vec<T> & x, const faust_vec<T> & y);
 
 // C = alpha *op(A)*op(B) + beta * C;
 // op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'
 // op(B) = B si typeB='N', op(B) = transpose(B) si typeB='T'
 // l'objet C doit etre different de A et B
template<typename T>
void gemm(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB);

template<typename T>
void gemv(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);

template<typename T> 
T power_iteration(const faust_mat<T> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag);

 // C = op(A) * B if typeMult = 'R'
 // C = B * op(A) sinon
 // op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'


// non-member operators declarations
template<typename T>
faust_vec<T> operator*(const faust_mat<T>& M, const faust_vec<T>& v);

#include "LinAlgebra.hpp"

#endif