/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef LINALGEBRA_H
#define LINALGEBRA_H

#include "faust_constant.h"
#include "faust_BlasHandle.h"
#include "faust_LinearOperator.h"
#include <complex>

namespace Faust
{

	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class MatGeneric;

	template <typename FPP>
		void setOp(const Faust::MatDense<FPP,Cpu>& A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);


	//modif AL AL
	//! \fn add
	//! \brief (*this) = (*this) + A
	template<typename FPP>
		void add(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C);


	//! \fn Faust::spgemm
	//! \brief performs Sparse matrices multiplication
	template<typename FPP>
		void spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);


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
		void gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB, Faust::BlasHandle<Cpu> const blas_handle)
		{Faust::gemm(A,B,C,alpha,beta,typeA,typeB);}


	//////////FONCTION Faust::MatDense<FPP,Cpu> - Faust::MatDense<FPP,Cpu> ////////////////////
	// C = A * B;
	//l'objet C doit etre different de A et B (but gemm() manages it with two copies)
	//! \fn multiply
	//! \brief Multiplication C = A * B
	//! \warning Object C must be different of A and B.
	template<typename FPP>
		void multiply(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C)
		{Faust::gemm(A,B,C,(FPP) 1.0,(FPP)0.0,'N','N');}
	template<typename FPP>
		void multiply(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C, Faust::BlasHandle<Cpu> const handle)
		{Faust::multiply(A,B,C);}

	// modif AL AL
	//modif AL AL
	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP>
		void multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
	// modif AL AL

	template<typename FPP>
		void gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);


	/*! \fn Faust::gemv
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
		void gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA,Faust::BlasHandle<Cpu> & handle)
		{Faust::gemv(A,x,y,alpha,beta,typeA);}



	/*!
	 *  \fn FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag)
	 *  \brief compute the biggest eigenvalue of A using power_iteration algorithm
	 * 	\param [in] A : input matrix must be semidefinite positive
	 *	\param [in] nbr_iter_max : maximum number of iteration
	 *	\param [in] threshold : threshold until convergence
	 *	\param [in,out] flag : convergence flag
	 *	\return	the the biggest eigenvalue of A
	 */
	//    template<typename FPP, typename FPP2 = double>
	//    FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP2 threshold,faust_int & flag);
	template<typename FPP, typename FPP2 = double>
		FPP power_iteration(const Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP2 threshold,faust_int & flag,Faust::BlasHandle<Cpu> & handle)
		{Faust::power_iteration(A,nbr_iter_max,threshold,flag);}
	    template<typename FPP, typename FPP2 = double>
		FPP power_iteration(const Faust::LinearOperator<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP2 threshold,int & flag);


	//! surcharge d'operateur * pour multiplier des matrices et des vecteurs
	template<typename FPP>
		Faust::Vect<FPP,Cpu> operator*(const Faust::MatDense<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v);
	template<typename FPP>
		Faust::Vect<FPP,Cpu> operator*(const Faust::MatSparse<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v);
	template<typename FPP>
		Faust::Vect<FPP,Cpu> operator*(const Faust::Transform<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v);


	//
	template<typename FPP>
		FPP fabs(FPP c);

	template<typename FPP>
		FPP fabs(std::complex<FPP> c);
	template<typename FPP, Device DEVICE>
		void multiply_order_opt_all_ends(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));
	template<typename FPP, Device DEVICE>
		void multiply_order_opt_all_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));
	template<typename FPP, Device DEVICE>
		void multiply_order_opt_first_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));
	template<typename FPP, Device DEVICE>
		void multiply_order_opt(const int mode, std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha=1.0, FPP beta_out=.0, std::vector<char> transconj_flags = std::vector<char>({'N'}));

}


#include "faust_linear_algebra.hpp"

#endif
