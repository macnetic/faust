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
#include <complex>
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
//			std::cout << "gemv identity opt." << std::endl;
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
	while((Faust::fabs(lambda_old-lambda)> Faust::fabs(threshold) || Faust::fabs(lambda) <= Faust::fabs(threshold)) && i<nbr_iter_max)
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
	//rewrite fabs even if c++11 standard provides it
	//because clang++ from macOS failed to provide fabs(complex)
	//(gcc-g++ provides it! but we compile with clang on macOS)
	return sqrt(norm(c));
}

	template<typename FPP>
void conjugate(std::complex<FPP>* elts, faust_unsigned_int n)
{
	for(faust_unsigned_int i=0; i< n; i++)
		elts[i] = std::complex<FPP>(elts[i].real(), - elts[i].imag());
}

	template<typename FPP>
void conjugate(FPP* elts, faust_unsigned_int n)
{
	//nothing to do for real numbers
}
#endif
