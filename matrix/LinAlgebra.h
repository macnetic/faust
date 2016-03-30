#ifndef LINALGEBRA_H
#define LINALGEBRA_H

#include "faust_constant.h"

template<typename T> class faust_vec;
template<typename T> class faust_mat;
template<typename T> class faust_spmat;

//////////FONCTION faust_mat<T> - faust_mat<T> ////////////////////
// C = A * B;
//l'objet C doit etre different de A et B
/*! \brief Multiplication C = A * B
 * Object C must be different of A and B.
 * */
template<typename T>
void multiply(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);


/*! \brief Gemm — Performs a matrix calculation of the form a(AB) +bC, with three matrices and two real constants you supply. <br>
* C = alpha *op(A)*op(B) + beta * C; <br>
* op(A) = A if typeA='N', op(A) = transpose(A) if typeA='T'<br>
* op(B) = B if typeB='N', op(B) = transpose(B) if typeB='T'<br>
* The Object C must be different of A and B. <br>
* \param A : Dense Matrix
* \param B : Dense Matrix
* \param C : Dense Matrix
* \param alpha : Template coefficient
* \param beta : Template coefficient
*/
template<typename T>
void gemm(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB);

/*! \brief Gemv — Performs a matrix vector multiplies of the form ????, with one matrice and two vectors and two real constants you supply. <br>
* C = alpha *op(A) + beta * C; <br>
* op(A) = A si typeA='N', op(A) = transpose(A) si typeA='T'<br>
* op(B) = B si typeB='N', op(B) = transpose(B) si typeB='T'<br>
* The Object C must be different of A and B. <br>
* \param A : Dense Matrix
* \param x : vector
* \param y : vector
* \param alpha : Template coefficient
* \param beta : Template coefficient
* \param typeA : char
*/
template<typename T>
void gemv(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);

/*!
*  \fn T power_iteration(const faust_mat<T> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag)
*  \brief compute the biggest eigenvalue of A using power_iteration algorithm
* 	\param[in] A : input matrix must be semidefinite positive
*	\param[in] nbr_iter_max : maximum number of iteration
*	\param[in] threshold : threshold until convergence
*	\param[in,out] flag : convergence flag
*	\return	the the biggest eigenvalue of A
*/
template<typename T>
T power_iteration(const faust_mat<T> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag);

template<typename T>
faust_vec<T> operator*(const faust_mat<T>& M, const faust_vec<T>& v);

#include "LinAlgebra.hpp"

#endif
