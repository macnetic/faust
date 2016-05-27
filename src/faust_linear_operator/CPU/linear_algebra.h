#ifndef LINALGEBRA_H
#define LINALGEBRA_H

#include "faust_constant.h"
#include "BlasHandleCPU.h"


template<typename FPP,Device DEVICE> class Vect;
template<typename FPP,Device DEVICE> class MatDense;
// modif AL AL



template <typename FPP>
void setOp(const Faust::MatDense<FPP,Cpu>& A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);



/*! \fn gemm
*   \brief Performs a matrix calculation of the form a(AB) +bC, with three matrices and two real constants you supply. <br>
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
template<typename FPP>
void gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
template<typename FPP>
void gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB,BlasHandle<Cpu> const blas_handle)
{gemm(A,B,C,alpha,beta,typeA,typeB);}


//////////FONCTION Faust::MatDense<FPP,Cpu> - Faust::MatDense<FPP,Cpu> ////////////////////
// C = A * B;
//l'objet C doit etre different de A et B
//! \fn multiply
//! \brief Multiplication C = A * B
//! \warning Object C must be different of A and B.
template<typename FPP>
void multiply(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C)
{gemm(A,B,C,(FPP) 1.0,(FPP)0.0,'N','N');}
template<typename FPP>
void multiply(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,BlasHandle<Cpu> const handle)
{multiply(A,B,C);}


template<typename FPP>
void gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);


/*! \fn gemv
*   \brief Performs a matrix vector multiplies of the form ????, with one matrice and two vectors and two real constants you supply. <br>
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
template<typename FPP>
void gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
template<typename FPP>
void gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA,BlasHandle<Cpu> & handle)
{gemv(A,x,y,alpha,beta,typeA);}



/*!
*  \fn FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag)
*  \brief compute the biggest eigenvalue of A using power_iteration algorithm
* 	\param [in] A : input matrix must be semidefinite positive
*	\param [in] nbr_iter_max : maximum number of iteration
*	\param [in] threshold : threshold until convergence
*	\param [in,out] flag : convergence flag
*	\return	the the biggest eigenvalue of A
*/
template<typename FPP>
FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP threshold,faust_int & flag);
template<typename FPP>
FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP threshold,faust_int & flag,BlasHandle<Cpu> & handle)
{power_iteration(A,nbr_iter_max,threshold,flag);}


//! surcharge d'operateur * pour multiplier des matrices et des vecteurs
template<typename FPP>
Faust::Vect<FPP,Cpu> operator*(const Faust::MatDense<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v);


#include "linear_algebra.hpp"

#endif
