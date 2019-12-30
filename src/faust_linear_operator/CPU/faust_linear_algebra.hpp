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
#ifndef __FAUST_LINALGEBRA_HPP
#define __FAUST_LINALGEBRA_HPP

#include <iostream>
#include <stdexcept>
#include "faust_LinearOperator.h"
#include "faust_MatGeneric.h"
#include "faust_MatDense.h"
#include "faust_constant.h"
#include "faust_Vect.h"

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif


#ifdef __GEMM_WITH_OPENBLAS__
	#include "faust_cblas_algebra.h"
#endif



	
template<typename FPP>
void Faust::spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB)
{
	//TODO: refactoring should be done to avoid repeating similar block of code for different cases (typeA,typeB,alpha,beta)
//#ifdef __COMPILE_TIMERS__
//	A.t_gemm.start();
//#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if (((&(C.mat)) == (&(B.mat))))
	{
		handleError("linear_algebra", " Faust::spgemm : C is the same object as A or B");
	}

	if (typeA == 'T' || typeA == 'H')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}


	if (typeB == 'T' || typeA == 'H')
	{
		nbRowOpB = B.getNbCol();
		nbColOpB = B.getNbRow();
	}else
	{
		nbRowOpB = B.getNbRow();
		nbColOpB = B.getNbCol();
	}


	if (nbColOpA != nbRowOpB)
	{
		handleError("linear_algebra", "Faust::spgemm : dimension conflict  between matrix op(A) and matrix op(B)");

	}


	if ( (beta!= FPP(0.0))  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		//handleError("Linalgebra : gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError("linear_algebra", "Faust::spgemm : invalid dimension for output matrix C");
	}

        C.resize(nbRowOpA,nbColOpB);





	if (beta == FPP(0.0))
	{

		if(B.isZeros)
		{

			FPP *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(FPP) * C.dim1*C.dim2);
			C.isZeros = true;
			C.set_id(false);
			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}

#define M_times_alpha_into_C(alpha, M, typeM) \
		{ \
			C=M; \
			if(typeM == 'T') \
			C.transpose(); \
			else if(typeM == 'H') \
			{ \
				C.transpose(); \
				C.conjugate(); \
			} \
			if(alpha!=FPP(1.0)) \
			C*= alpha; \
			return; \
		}


		if(B.is_id())
			M_times_alpha_into_C(alpha, A, typeA); //return here

		if(A.is_id())
			M_times_alpha_into_C(alpha, B, typeB); //return here

		if (typeA == 'N')
		{
			if (typeB == 'N')
				C.mat.noalias() = alpha * A.mat * B.mat;
			else if(typeB == 'T')
				C.mat.noalias() = alpha * A.mat * B.mat.transpose();
			else // typeB == 'H'
				C.mat.noalias() = alpha * A.mat * B.mat.adjoint();
		}else if(typeA == 'T')
		{
			if (typeB == 'N')
				C.mat.noalias() = alpha * A.mat.transpose() * B.mat;
			else if(typeB == 'T')
				C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose();
			else // typeB == 'H'
				C.mat.noalias() = alpha * A.mat.transpose() * B.mat.adjoint();
		} else // typeA == 'H'
		{
			if (typeB == 'N')
				C.mat.noalias() = alpha * A.mat.adjoint() * B.mat;
			else if(typeB == 'T')
				C.mat.noalias() = alpha * A.mat.adjoint() * B.mat.transpose();
			else // typeB == 'H'
				C.mat.noalias() = alpha * A.mat.adjoint() * B.mat.adjoint();

		}


	}else //beta != 0
	{
		if(B.isZeros)
		{
			C *= beta;
			C.isZeros = false;
			C.set_id(false);
			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}

#define M_times_alpha_plus_beta_times_C_into_C(M, typeM, alpha, beta, C)\
		{ \
			C *= beta; \
			if(typeM == 'N' && alpha == FPP(1.0)) \
			{ \
				C += M; \
				C.isZeros = false; \
				C.set_id(false); \
			} \
			Faust::MatDense<FPP,Cpu> M_tmp(M); \
			if(typeM == 'T') \
			M_tmp.transpose(); \
			else if(typeM == 'H') \
			{ \
				M_tmp.conjugate(false); \
				M_tmp.transpose(); \
			} \
			if(alpha != FPP(1.0)) \
			M_tmp *= alpha; \
			C += M_tmp; \
			return; \
		}

		if(B.is_id())
			M_times_alpha_plus_beta_times_C_into_C(A, typeA, alpha, beta, C); // return here

		if(A.is_id())
			M_times_alpha_plus_beta_times_C_into_C(B, typeB, alpha, beta, C); //return here

		if (typeA == 'N')
		{
			if (typeB == 'N')
				C.mat = alpha * A.mat * B.mat + beta * C.mat;
			else if(typeB == 'T')
				C.mat = alpha * A.mat * B.mat.transpose() + beta * C.mat;
			else //typeB == 'H'
				C.mat = alpha * A.mat * B.mat.adjoint() + beta * C.mat;
		}else if(typeA == 'T')
		{
			if (typeB == 'N')
				C.mat = alpha * A.mat.transpose() * B.mat + beta * C.mat ;
			else if(typeB == 'T')
				C.mat = alpha * A.mat.transpose() * B.mat.transpose() + beta * C.mat;
			else //typeB 'H'
				C.mat = alpha * A.mat.transpose() * B.mat.adjoint() + beta * C.mat;
		}
		else //typeA == 'H'
		{
			if (typeB == 'N')
				C.mat = alpha * A.mat.adjoint() * B.mat + beta * C.mat ;
			else if(typeB == 'T')
				C.mat = alpha * A.mat.adjoint() * B.mat.transpose() + beta * C.mat;
			else //typeB 'H'
				C.mat = alpha * A.mat.adjoint() * B.mat.adjoint() + beta * C.mat;
		}

	}
	C.isZeros = false;
	C.set_id(false);
//#ifdef __COMPILE_TIMERS__
//A.t_gemm.stop();
//#endif
}

template<typename FPP>
void Faust::gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA)
{
	faust_unsigned_int nbRowOpA,nbColOpA;
	const Faust::Vect<FPP,Cpu>* px;
	if  ((&x) == (&y))
		px = new Faust::Vect<FPP,Cpu>(x);
	else
		px = &x;


	if (typeA == 'T' || typeA == 'H')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}

	if   (nbColOpA != px->size() )
	{
		//handleError("Linalgebra : Faust::gemv : nbCol of op(A) = %d while dim of x = %d",nbColOpA,px->getDim());
		handleError("linear_algebra", "Faust::gemv : dimension conflict  between matrix op(A) and input std::vector x");
	}

	if ( (beta!=FPP(0.0))  &&  (y.size() != nbRowOpA))
	{
		handleError("linear_algebra", "Faust::gemv : dimension conflict  between matrix op(A) and output std::vector y");
	}

	y.resize(nbRowOpA);


	#ifdef __GEMM_WITH_OPENBLAS__
		CBLAS_TRANSPOSE transA,transB;
		if (typeA=='T' || typeB == 'H')
			transA = CblasTrans;
		else
			transA = CblasNoTrans;
	#endif

	#ifndef __GEMM_WITH_OPENBLAS__
	if (beta == FPP(0.0))
	{
		if(A.is_id() && alpha == FPP(1.0))
		{
			std::cout << "gemv identity opt." << std::endl;
			memcpy(y.getData(), px->getData(), sizeof(FPP)*nbRowOpA);
		}
		else if (typeA == 'N')
		{
			y.vec.noalias() = alpha * A.mat * px->vec;
		}else if (typeA == 'T')
		{

			y.vec.noalias() = alpha * A.mat.transpose() * px->vec;
		}else // typeA == 'H'
		{
			y.vec.noalias() = alpha * A.mat.adjoint() * px->vec;
		}
	}else
	{
		if (typeA == 'N')
		{
			y.vec = alpha * A.mat * px->vec + beta * y.vec;
		}else if(typeA == 'T')
		{
			y.vec = alpha * A.mat.transpose() * px->vec + beta * y.vec;
		}
		else // typeA == 'H'
		{
			y.vec = alpha * A.mat.adjoint() * px->vec + beta * y.vec;
		}
	}
	#else

		Faust::cblas_gemv<FPP>(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);
	#endif

	if  ((&x) == (&y))
		delete px;
	px=NULL;

}

template<typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
{
	if ( (&C == &A) || (&C == &B) )
	{
		Faust::MatDense<FPP,Cpu> Cbis(C);
		Faust::gemm_core(A,B,Cbis,alpha,beta,typeA,typeB);
		C=Cbis;
	}else
	{Faust::gemm_core(A,B,C,alpha,beta,typeA,typeB);}
}


//WARNING matrix C must be a different object from A and B
template<typename FPP>
void Faust::gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
{

#ifdef __COMPILE_TIMERS__
	A.t_gemm.start();
#endif

#ifdef __GEMM_WITH_OPENBLAS__
	if(typeA == 'H' || typeB == 'H')
		handleError("linear_algebra", " gemm: complex adjoint matrix is not yet handled with BLAS.");
#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{

		handleError("linear_algebra", " gemm : C is the same object as A or B");
	}

	if (typeA == 'T' || typeA == 'H')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}


	if (typeB == 'T' || typeB == 'H')
	{
		nbRowOpB = B.getNbCol();
		nbColOpB = B.getNbRow();
	}else
	{
		nbRowOpB = B.getNbRow();
		nbColOpB = B.getNbCol();
	}



	if (nbColOpA != nbRowOpB)
	{
		handleError("linear_algebra", "Faust::gemm : dimension conflict  between matrix op(A) and matrix op(B)");

	}


	if ( (beta!=FPP(0.0))  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		//handleError("Linalgebra : Faust::gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError("linear_algebra", "Faust::gemm : invalid dimension for output matrix C");
	}

        C.resize(nbRowOpA,nbColOpB);


	#ifdef __GEMM_WITH_OPENBLAS__
		CBLAS_TRANSPOSE transA,transB;
		if (typeA=='T')
			transA = CblasTrans;
		else
			transA = CblasNoTrans;
		if (typeB=='T')
			transB = CblasTrans;
		else
			transB = CblasNoTrans;
	#endif



	if (beta == FPP(0.0))
	{

		if(A.isZeros || B.isZeros)
		{

			FPP *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(FPP) * C.dim1*C.dim2);
			C.isZeros = true;
			C.set_id(false);
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}

#define N_times_alpha_into_C(N, typeN, alpha, C) \
		{ \
			C = N; \
			if(typeN == 'T') \
			C.transpose(); \
			else if(typeN == 'H') \
			{ \
				C.conjugate(false); \
				C.transpose(); \
			} \
			if(alpha!=FPP(1.0)) \
			C *= alpha; \
			return; \
		}

		if(A.is_id())
			N_times_alpha_into_C(B, typeB, alpha, C);

		if(B.is_id())
			N_times_alpha_into_C(A, typeA, alpha, C);

		#ifndef __GEMM_WITH_OPENBLAS__
		// std::cout<<" A normal Faust::gemm"<<std::endl;
			if (typeA == 'N')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat * B.mat;
				else if(typeB == 'T')
					C.mat.noalias() = alpha * A.mat * B.mat.transpose();
				else // typeB == 'H' //TODO: check validity of typeA and typeB
					C.mat.noalias() = alpha * A.mat * B.mat.transpose().conjugate();
			}else if(typeA == 'T')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat;
				else if(typeB == 'T')
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose();
				else // typeB == 'H' //TODO: check validity of typeA and typeB
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose().conjugate();
			}
			else //if(typeA == 'H')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat.transpose().conjugate() * B.mat;
				else if(typeB == 'T')
					C.mat.noalias() = alpha * A.mat.transpose().conjugate() * B.mat.transpose();
				else // typeB == 'H' //TODO: check validity of typeA and typeB
					C.mat.noalias() = alpha * A.mat.transpose().conjugate() * B.mat.transpose().conjugate();
			}
		#else
			 FPP beta = FPP(0.0);
			 Faust::cblas_gemm<FPP>(CblasColMajor, transA, transB, (int) C.dim1, (int)  C.dim2, (int) nbColOpA, (FPP) alpha, (FPP*) A.getData(), (int) A.dim1, (FPP*) B.getData(), (int) B.dim1,(FPP) beta, (FPP*) C.getData(),(int) C.dim1);

		#endif


	}else
	{
		if(A.isZeros || B.isZeros)
		{
			C *= beta;
			C.isZeros = false;
			C.set_id(false);
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}

		if(A.is_id())
		{
			C *= beta;
			if(typeB == 'N' && alpha == FPP(1.0))
			{
				C += B;
				#ifdef __COMPILE_TIMERS__
					A.t_gemm.stop();
				#endif
				return;
			}
			Faust::MatDense<FPP,Cpu> B_tmp(B);
			if(typeB == 'T')
				B_tmp.transpose();
			else if(typeB == 'H')
			{
				B_tmp.conjugate(false);
				B_tmp.transpose();
			}
			if(alpha != FPP(1.0))
				B_tmp *= alpha;
			C += B_tmp;
			C.isZeros = false;
			C.set_id(false);
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
		if(B.is_id())
		{
			C *= beta;
			if(typeA == 'N' && alpha == FPP(1.0))
			{
				C += A;
				C.isZeros = false;
				C.set_id(false);
				#ifdef __COMPILE_TIMERS__
					A.t_gemm.stop();
				#endif
				return;
			}
			Faust::MatDense<FPP,Cpu> A_tmp(A);
			if(typeA == 'T')
				A_tmp.transpose();
			else if(typeA == 'H')
			{
				A_tmp.conjugate(false);
				A_tmp.transpose();
			}
			if(alpha != FPP(1.0))
				A_tmp *= alpha;
			C += A_tmp;
			C.isZeros = false;
			C.set_id(false);

			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}


		#ifndef __GEMM_WITH_OPENBLAS__
			if (typeA == 'N')
			{
				if (typeB == 'N')
						C.mat = alpha * A.mat * B.mat + beta * C.mat;
				else
					C.mat = alpha * A.mat * B.mat.transpose() + beta * C.mat;
			}else
			{
				if (typeB == 'N')
					C.mat = alpha * A.mat.transpose() * B.mat + beta * C.mat ;
				else
					C.mat = alpha * A.mat.transpose() * B.mat.transpose() + beta * C.mat;
			}
		#else
			Faust::cblas_gemm<FPP>(CblasColMajor, transA, transB, (int) C.dim1,(int)  C.dim2,(int)  nbColOpA,(FPP) alpha,(FPP*)  A.getData(), (int) A.dim1,(FPP*) B.getData(),(int) B.dim1, (FPP) beta, (FPP*) C.getData(), (int)C.dim1);



		#endif

	}
	C.isZeros = false;
	C.set_id(false);
#ifdef __COMPILE_TIMERS__
A.t_gemm.stop();
#endif
}



template<typename FPP>
void Faust::add(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C)
{
#ifdef __COMPILE_TIMERS__
A.t_add_ext.start();
#endif
	if ((A.getNbCol() != B.getNbCol()) || (A.getNbRow() != B.getNbRow()) || (A.getNbRow() != C.getNbRow()) || (A.getNbCol() != C.getNbCol()))
	{
		handleError("linear_algebra"," add : matrix dimension not equal");
	}else
	{
		C.mat = A.mat + B.mat;
	}
	C.isZeros = false;
	C.set_id(false);
#ifdef __COMPILE_TIMERS__
A.t_add_ext.stop();
#endif
}



// compute the biggest eigenvalue of A, A must be semi-definite positive
template<typename FPP, typename FPP2>
FPP Faust::power_iteration(const  Faust::LinearOperator<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP2 threshold, int & flag)
{
	#ifdef __COMPILE_TIMERS__
		A.t_power_iteration.start();
	#endif


	const int nb_col = A.getNbCol();
	int i = 0;
	flag = 0;

	if (nbr_iter_max <= 0)
    {
        #ifdef __COMPILE_TIMERS__
            A.t_power_iteration.stop();
        #endif
   		handleError("linear_algebra "," power_iteration :  nbr_iter_max <= 0");
    }
	if (nb_col != A.getNbRow())
	{
	    #ifdef __COMPILE_TIMERS__
            A.t_power_iteration.stop();
        #endif
        handleError("linear_algebra "," power_iteration : Faust::Transform<FPP,Cpu> 1 must be a squared matrix");
	}
	Faust::Vect<FPP,Cpu> xk(nb_col);
	xk.setOnes();
	Faust::Vect<FPP,Cpu> xk_norm(nb_col);
	FPP lambda_old=1.0;
   	FPP lambda = 0.0;
   	FPP alpha = 1.0;
	FPP beta = 0.0;
	while(Faust::fabs(lambda_old-lambda)> Faust::fabs(threshold) && i<nbr_iter_max)
	{
		i++;
      		lambda_old = lambda;
      		xk_norm = xk;
      		xk_norm.normalize();
      		xk = A.multiply(xk_norm);
      		lambda = xk_norm.dot(xk);
     		//std::cout << "i = " << i << " ; lambda=" << lambda << std::endl;
   	}
   	flag = (i<nbr_iter_max)?i:-1;

    #ifdef __COMPILE_TIMERS__
        A.t_power_iteration.stop();
	#endif

   	return lambda;



}



// non-member operators definitions
template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::operator*(const Faust::MatDense<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v)
{
	Faust::Vect<FPP,Cpu> vec(v);
	vec.multiplyLeft(M);
	return vec;
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::operator*(const Faust::MatSparse<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v)
{
	Faust::Vect<FPP,Cpu> vec(v);
	vec.multiplyLeft(M);
	return vec;
}




template<typename FPP>
FPP Faust::fabs(FPP f)
{
	return std::fabs(f);
}

template<typename FPP>
FPP Faust::fabs(std::complex<FPP> c)
{
	//rewrite fabs even if c++11 standard provide it
	//because clang++ from macOS failed to provide fabs(complex) 
	//(gcc-g++ provides it! but we compile with clang on macOS)
	return sqrt(norm(c));
}

//TODO: allow nflags == 1 and all factors using with same flag
/**
 *	\brief Multiply all the matrices together (as a product of n factors) with an optimization based on the associativity of the matrix product ; following an order that minimizes the number of scalar operations.
 *	\note: the std::vector facts is altered after the call! Don't reuse it.
 */
template<typename FPP, Device DEVICE>
void Faust::multiply_order_opt_all_ends(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	Faust::MatDense<FPP, DEVICE> tmpr, tmpl;
	int nfacts = facts.size();
	int ri = 0, li = nfacts-1;
	Faust::MatDense<FPP,DEVICE> *R1, *R2, *L1, *L2;
	faust_unsigned_int R1nr, R1nc, L1nr, L1nc;
	while(li-ri > 1)
	{
		R1 = facts[ri];
		R2 = facts[ri+1];
		L1 = facts[li-1];
		L2 = facts[li];
		R1nr = R1->getNbRow();
		R1nc = R1->getNbCol();
		L1nr = L1->getNbRow();
		L1nc = L1->getNbCol();
		if(R1nr * R1nc * R2->getNbCol() < L1nr * L1nc * L2->getNbCol())
		{
			gemm(*R1, *R2, tmpr, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>ri?ri:0], transconj_flags[transconj_flags.size()>ri+1?ri+1:0]);
			ri++;
			facts[ri] = &tmpr;
			if(transconj_flags.size() > ri)
				transconj_flags[ri] = 'N';
		}
		else
		{
			gemm(*L1, *L2, tmpl, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>li-1?li-1:0], transconj_flags[transconj_flags.size()>li?li:0]);
			li--;
			facts[li] = &tmpl;
			if(transconj_flags.size() > li)
				transconj_flags[li] = 'N';
		}
	}
	// last mul
	gemm(*facts[ri], *facts[li], out, alpha, beta_out, ri==0?transconj_flags[0]:'N', li==nfacts-1&&transconj_flags.size()>li?transconj_flags[li]:'N');
	facts.erase(facts.begin(), facts.end());
}

template<typename FPP, Device DEVICE>
void Faust::multiply_order_opt_all_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	std::vector<Faust::MatDense<FPP,DEVICE>*> tmp_facts; //temporary product results
	Faust::MatDense<FPP, DEVICE>* tmp;
	int nfacts = facts.size();
	Faust::MatDense<FPP,DEVICE> *Si, *Sj;
	std::vector<int> complexity(nfacts-1);
	for(int i = 0; i <nfacts-1; i++)
	{
		Si = facts[i];
		Sj = facts[i+1];
		complexity[i] = Si->getNbRow() * Si->getNbCol() * Sj->getNbCol();
	}
	int idx; // marks the factor to update with a product of contiguous factors
	bool multiplying_tmp_factor = false; // allows to avoid to allocate uselessly a tmp factor if Si or Sj are already tmp factors
	while(facts.size() > 2)
	{
		// find the least complex product facts[idx]*facts[idx+1]
		idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
		Si = facts[idx];
		Sj = facts[idx+1];
		for(auto Tit = tmp_facts.begin(); Tit != tmp_facts.end(); Tit++)
		{
			if(Sj == *Tit)
			{// Sj is original fact
				multiplying_tmp_factor = true;
				tmp = Sj;
				break;
			}
			else if(Si == *Tit)
			{
				multiplying_tmp_factor = true;
				tmp = Si;
				break;
			}
		}
		if(! multiplying_tmp_factor)
		{
			tmp = new Faust::MatDense<FPP, DEVICE>();
			tmp_facts.push_back(tmp);
		}
		//else no need to instantiate a new tmp, erasing Sj which is a tmp
		gemm(*Si, *Sj, *tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>idx?idx:0], transconj_flags[transconj_flags.size()>idx+1?idx+1:0]);
		facts.erase(facts.begin()+idx+1);
		complexity.erase(complexity.begin()+idx); //complexity size == facts size - 1
		facts[idx] = tmp;
		if(transconj_flags.size() > idx)
			transconj_flags[idx] = 'N';
		// update complexity around the new factor
		if(facts.size() > 2)
		{
			if(idx > 0)
				complexity[idx-1] = facts[idx-1]->getNbRow() * facts[idx-1]->getNbCol() * facts[idx]->getNbCol();
			if(idx < facts.size()-1)
				complexity[idx] = facts[idx]->getNbRow() * facts[idx]->getNbCol() * facts[idx+1]->getNbCol();
		}
		multiplying_tmp_factor = false;
	}
	// last mul
	gemm(*facts[0], *facts[1], out, alpha, beta_out, transconj_flags[0], transconj_flags.size()>1?transconj_flags[1]:'N');
	facts.erase(facts.begin(), facts.end());
	// delete all tmp facts
	for(auto Tit = tmp_facts.begin(); Tit != tmp_facts.end(); Tit++)
	{
		delete *Tit;
	}
}

template<typename FPP, Device DEVICE>
void Faust::multiply_order_opt_first_best(std::vector<Faust::MatDense<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	int nfacts = facts.size();
	Faust::MatDense<FPP,DEVICE> *Si, *Sj, tmp;
	if(nfacts == 1)
	{
		tmp = *facts[0];
		if(transconj_flags[0] != 'N')
		{
			tmp.conjugate(false);
			tmp.transpose();
		}
		tmp *= alpha;
		if(beta_out != FPP(0))
		{
			out *= beta_out;
			out += tmp;
		}
		else
			out = tmp;

	}
	if(nfacts <= 2)
	{
		gemm(*facts[0], *facts[1], out, alpha, beta_out, transconj_flags[0], transconj_flags[transconj_flags.size()>1?1:0]);
		return;
	}
	std::vector<int> complexity(nfacts-1);
	int min_cplx = std::numeric_limits<int>::max();
	int i, idx, lasti; // idx marks the factor to update with a product of contiguous factors
	for(int i = 0; i <nfacts-1; i++)
	{
		Si = facts[i];
		Sj = facts[i+1];
		complexity[i] = Si->getNbRow() * Si->getNbCol() * Sj->getNbCol();
		if(complexity[i] < min_cplx)
		{
			min_cplx = complexity[i];
			idx = i;
		}
	}
	// find the least complex product facts[idx]*facts[idx+1]
//	idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
//	std::cout << "idx: " << idx << std::endl;
	Si = facts[idx];
	Sj = facts[idx+1];
	gemm(*Si, *Sj, tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>idx?idx:0], transconj_flags[transconj_flags.size()>idx+1?idx+1:0]);
	lasti = idx+2<facts.size()?0:1;
	i = idx-1;
	while(i >= lasti)
	{
		gemm(*facts[i], tmp, tmp, (FPP)1.0, (FPP)0.0, transconj_flags[transconj_flags.size()>i?i:0], 'N');
		i--;
	}
	if(lasti == 0)
	{
		// all left factors multiplied
		// now multiply on the right
		i = idx+2;
		while(i < facts.size()-1)
		{
			gemm(tmp, *facts[i], tmp, (FPP)1.0, (FPP)0, 'N', transconj_flags[transconj_flags.size()>i?i:0]);
			i++;
		}
		//last mul
		gemm(tmp, *facts[i], out, alpha, beta_out, 'N', transconj_flags[transconj_flags.size()>i?i:0]);
	}
	else
		gemm(*facts[i], tmp, out, alpha, beta_out, transconj_flags[transconj_flags.size()>i?i:0], 'N');
}

template<typename FPP, Device DEVICE>
void Faust::multiply_order_opt(const int mode, std::vector<Faust::MatGeneric<FPP,DEVICE>*>& facts, Faust::MatDense<FPP,DEVICE>& out, FPP alpha/* =1.0*/, FPP beta_out/*=.0*/, std::vector<char> transconj_flags /* = {'N'}*/)
{
	std::vector<Faust::MatDense<FPP,DEVICE>*> dfacts;
	std::vector<Faust::MatDense<FPP,DEVICE>*> dfacts_to_del;
	Faust::MatDense<FPP,DEVICE> * df = nullptr;
	for(auto f: facts)
	{
		if(!(df = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(f)))
		{
			Faust::MatSparse<FPP,DEVICE> *sf = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(f);
			df = new Faust::MatDense<FPP,DEVICE>(*sf);
			dfacts_to_del.push_back(df);
		}
		dfacts.push_back(df);
	}
	switch(mode)
	{
		case 1:
			Faust::multiply_order_opt_all_ends(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		case 2:
			Faust::multiply_order_opt_first_best(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		case 3:
			Faust::multiply_order_opt_all_best(dfacts, out, alpha, beta_out, transconj_flags);
			break;
		default:
			throw std::out_of_range("optimization asked not known");
	}
	// free dense copies
	for(auto df: dfacts_to_del)
		delete df;
}
#endif
