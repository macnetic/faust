#ifndef __FAUST_LINALGEBRA_HPP
#define __FAUST_LINALGEBRA_HPP

#include <iostream>

#include "faust_Vect.h"
#include "faust_MatDense.h"
#include "faust_constant.h"



#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif


#ifdef __GEMM_WITH_OPENBLAS__
	#include "faust_cblas_algebra.h"
#endif




template<typename FPP>
void spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB)
{
//#ifdef __COMPILE_TIMERS__
//	A.t_gemm.start();
//#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if (((&(C.mat)) == (&(B.mat))))
	{
		handleError("linear_algebra", " spgemm : C is the same object as A or B");
	}

	if (typeA == 'T')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}


	if (typeB == 'T')
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
		handleError("linear_algebra", "spgemm : dimension conflict  between matrix op(A) and matrix op(B)");

	}


	if ( (beta!=0)  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		//handleError("Linalgebra : gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError("linear_algebra", "spgemm : invalid dimension for output matrix C");
	}

        C.resize(nbRowOpA,nbColOpB);





	if (beta == 0.0)
	{

		if(B.isZeros)
		{

			FPP *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(FPP) * C.dim1*C.dim2);
			C.isZeros = true;
			C.isIdentity = false;
			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}
		if(B.isIdentity)
		{
			C=A;
			if(typeA == 'T')
				C.transpose();
			if(alpha!=1.0)
				C*= alpha;
			C.isZeros = false;
			C.isIdentity = false;
			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}

			if (typeA == 'N')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat * B.mat;
				else
					C.mat.noalias() = alpha * A.mat * B.mat.transpose();
			}else
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat;
				else
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose();
			}


	}else
	{
		if(B.isZeros)
		{
			C *= beta;
			C.isZeros = false;
			C.isIdentity = false;
			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}
		if(B.isIdentity)
		{
			C *= beta;
			if(typeA == 'N' && alpha == 1.0)
			{
				C += A;
				C.isZeros = false;
				C.isIdentity = false;
				//#ifdef __COMPILE_TIMERS__
				//	A.t_gemm.stop();
				//#endif
				return;
			}
			Faust::MatDense<FPP,Cpu> A_tmp(A);
			if(typeA == 'T')
				A_tmp.transpose();
			if(alpha != 1.0)
				A_tmp *= alpha;
			C += A_tmp;
			C.isZeros = false;
			C.isIdentity = false;

			//#ifdef __COMPILE_TIMERS__
			//	A.t_gemm.stop();
			//#endif
			return;
		}

			if (typeA == 'N')
			{
				if (typeB == 'N')
				{
						C.mat = alpha * A.mat * B.mat + beta * C.mat;
				}else
					C.mat = alpha * A.mat * B.mat.transpose() + beta * C.mat;
			}else
			{
				if (typeB == 'N')
					C.mat = alpha * A.mat.transpose() * B.mat + beta * C.mat ;
				else
					C.mat = alpha * A.mat.transpose() * B.mat.transpose() + beta * C.mat;
			}


	}
	C.isZeros = false;
	C.isIdentity = false;
//#ifdef __COMPILE_TIMERS__
//A.t_gemm.stop();
//#endif
}

template<typename FPP>
void gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA)
{
	faust_unsigned_int nbRowOpA,nbColOpA;
	const Faust::Vect<FPP,Cpu>* px;
	if  ((&x) == (&y))
		px = new Faust::Vect<FPP,Cpu>(x);
	else
		px = &x;


	if (typeA == 'T')
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
		//handleError("Linalgebra : gemv : nbCol of op(A) = %d while dim of x = %d",nbColOpA,px->getDim());
		handleError("linear_algebra", "gemv : dimension conflict  between matrix op(A) and input vector x");
	}

	if ( (beta!=0)  &&  (y.size() != nbRowOpA))
	{
		handleError("linear_algebra", "gemv : dimension conflict  between matrix op(A) and output vector y");
	}

	y.resize(nbRowOpA);


	#ifdef __GEMM_WITH_OPENBLAS__
		CBLAS_TRANSPOSE transA,transB;
		if (typeA=='T')
			transA = CblasTrans;
		else
			transA = CblasNoTrans;
	#endif

	#ifndef __GEMM_WITH_OPENBLAS__
	if (beta == 0.0)
	{
		if (typeA == 'N')
		{
			y.vec.noalias() = alpha * A.mat * px->vec;
		}else
		{

			y.vec.noalias() = alpha * A.mat.transpose() * px->vec;
		}
	}else
	{
		if (typeA == 'N')
		{
			y.vec = alpha * A.mat * px->vec + beta * y.vec;
		}else
		{
			y.vec = alpha * A.mat.transpose() * px->vec + beta * y.vec;
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
void gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
{
	if ( (&C == &A) || (&C == &B) )
	{
		Faust::MatDense<FPP,Cpu> Cbis(C);
		gemm_core(A,B,Cbis,alpha,beta,typeA,typeB);
		C=Cbis;
	}else
	{gemm_core(A,B,C,alpha,beta,typeA,typeB);}
}


//WARNING matrix C must be a different object from A and B
template<typename FPP>
void gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
{

#ifdef __COMPILE_TIMERS__
	A.t_gemm.start();
#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{

		handleError("linear_algebra", " gemm : C is the same object as A or B");
	}

	if (typeA == 'T')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}


	if (typeB == 'T')
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
		handleError("linear_algebra", "gemm : dimension conflict  between matrix op(A) and matrix op(B)");

	}


	if ( (beta!=0)  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		//handleError("Linalgebra : gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError("linear_algebra", "gemm : invalid dimension for output matrix C");
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



	if (beta == 0.0)
	{

		if(A.isZeros || B.isZeros)
		{

			FPP *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(FPP) * C.dim1*C.dim2);
			C.isZeros = true;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}

		if(A.isIdentity)
		{
			C=B;
			if(typeB == 'T')
				C.transpose();
			if(alpha!=1.0)
				C*= alpha;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
		if(B.isIdentity)
		{
			C=A;
			if(typeA == 'T')
				C.transpose();
			if(alpha!=1.0)
				C*= alpha;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}


		#ifndef __GEMM_WITH_OPENBLAS__
		// std::cout<<" A normal gemm"<<std::endl;
			if (typeA == 'N')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat * B.mat;
				else
					C.mat.noalias() = alpha * A.mat * B.mat.transpose();
			}else
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat;
				else
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose();
			}
		#else
			 FPP beta = 0.0;
			 Faust::cblas_gemm<FPP>(CblasColMajor, transA, transB, (int) C.dim1, (int)  C.dim2, (int) nbColOpA, (FPP) alpha, (FPP*) A.getData(), (int) A.dim1, (FPP*) B.getData(), (int) B.dim1,(FPP) beta, (FPP*) C.getData(),(int) C.dim1);

		#endif


	}else
	{
		if(A.isZeros || B.isZeros)
		{
			C *= beta;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}

		if(A.isIdentity)
		{
			C *= beta;
			if(typeB == 'N' && alpha == 1.0)
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
			if(alpha != 1.0)
				B_tmp *= alpha;
			C += B_tmp;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
		if(B.isIdentity)
		{
			C *= beta;
			if(typeA == 'N' && alpha == 1.0)
			{
				C += A;
				C.isZeros = false;
				C.isIdentity = false;
				#ifdef __COMPILE_TIMERS__
					A.t_gemm.stop();
				#endif
				return;
			}
			Faust::MatDense<FPP,Cpu> A_tmp(A);
			if(typeA == 'T')
				A_tmp.transpose();
			if(alpha != 1.0)
				A_tmp *= alpha;
			C += A_tmp;
			C.isZeros = false;
			C.isIdentity = false;

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
	C.isIdentity = false;
#ifdef __COMPILE_TIMERS__
A.t_gemm.stop();
#endif
}



template<typename FPP>
void add(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C)
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
	C.isIdentity = false;
#ifdef __COMPILE_TIMERS__
A.t_add_ext.stop();
#endif
}



// compute the biggest eigenvalue of A, A must be semi-definite positive
template<typename FPP>
FPP power_iteration(const  Faust::MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag)
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
	while(fabs(lambda_old-lambda)>threshold && i<nbr_iter_max)
	{
		i++;
      		lambda_old = lambda;
      		xk_norm = xk;
      		xk_norm.normalize();
      		xk = A*xk_norm;
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
Faust::Vect<FPP,Cpu> operator*(const Faust::MatDense<FPP,Cpu>& M, const Faust::Vect<FPP,Cpu>& v)
{
	Faust::Vect<FPP,Cpu> vec(v);
	vec.multiplyLeft(M);
	return vec;
}

#endif
