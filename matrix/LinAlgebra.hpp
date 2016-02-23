#ifndef __FAUST_LINALGEBRA_HPP
#define __FAUST_LINALGEBRA_HPP

#include <iostream>

#include "faust_vec.h"
#include "faust_mat.h"
#include "faust_constant.h"


	//////////FONCTION faust_mat<T> - faust_mat<T> ////////////////////

#ifdef __COMPILE_TIMERS__
	#include "faust_timer.h"
#endif


#ifdef __GEMM_WITH_OPENBLAS__
	#include "cblas_algebra.h"
#endif


/*template <typename T>
T dot(const faust_vec<T>& v1, const faust_vec<T>& v2)
{
   if(v1.size() != v2.size())
      handleError("LinAlgebra","dot : the two vectors don't have the same size");
      T result = v1.vec.dot(v2.vec);
      return result;
}*/





template<typename T>
 void multiply(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C)
{   
#ifdef __COMPILE_TIMERS__
A.t_multiply.start();
#endif

	if (A.getNbCol() != B.getNbRow())
	{
		//handleError("Linalgebra : multiply :  nbCol of A = %d while nbRow of B = %d",A.getNbCol(),B.getNbRow());
		handleError("LinAlgebra","multiply : invalid matrix dimensions");	
	}
	 
	if(A.isZeros || B.isZeros)
	{
		C.resize(A.dim1,B.dim2);
		T *const ptr_data_dst = C.getData();
		memset(ptr_data_dst, 0, sizeof(T) * C.dim1*C.dim2);
		C.isZeros = true;
		C.isIdentity = false;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}

	if(B.isIdentity)
	{
		C=A;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}

	if(A.isIdentity)
	{
		C=B;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}


	if (((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))))
	{
		handleError("LinAlgebra"," multiply : C is the same object as A or B");		
	}else
	{
		C.resize(A.getNbRow(),B.getNbCol());
		C.mat.noalias() = A.mat * B.mat;
	}
	C.isZeros = false;
	C.isIdentity = false;

#ifdef __COMPILE_TIMERS__
A.t_multiply.stop();
#endif
}



template<typename T>
void gemv(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA)
{
	faust_unsigned_int nbRowOpA,nbColOpA;
	const faust_vec<T>* px;
	if  ((&x) == (&y))
		px = new faust_vec<T>(x);
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
		handleError("LinAlgebra", "gemv : dimension conflict  between matrix op(A) and input vector x");
	}
	
	if ( (beta!=0)  &&  (y.size() != nbRowOpA))
	{
		handleError("LinAlgebra", "gemv : dimension conflict  between matrix op(A) and output vector y");
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

		cblas_gemv<T>(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);
	#endif
	
	
							
	if  ((&x) == (&y))
		delete px; 
	px=NULL;
	
}
	
	

template<typename T>	
void gemm(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB)
{

#ifdef __COMPILE_TIMERS__
A.t_gemm.start();
#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{
		handleError("LinAlgebra", "gemv : gemm : C is the same object as A or B");		
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
		handleError("LinAlgebra", "gemm : dimension conflict  between matrix op(A) and matrix op(B)");
		
	}







	if ( (beta!=0)  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{

		//handleError("Linalgebra : gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError("LinAlgebra", "gemm : invalid dimension for output matrix C");
	
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
			
			T *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(T) * C.dim1*C.dim2);
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
			 T beta = 0.0;	
			 cblas_gemm<T>(CblasColMajor, transA, transB, (int) C.dim1, (int)  C.dim2, (int) nbColOpA, (T) alpha, (T*) A.getData(), (int) A.dim1, (T*) B.getData(), (int) B.dim1,(T) beta, (T*) C.getData(),(int) C.dim1);

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
			faust_mat<T> B_tmp(B);
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
			faust_mat<T> A_tmp(A);
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
			cblas_gemm<T>(CblasColMajor, transA, transB, (int) C.dim1,(int)  C.dim2,(int)  nbColOpA,(T) alpha,(T*)  A.getData(), (int) A.dim1,(T*) B.getData(),(int) B.dim1, (T) beta, (T*) C.getData(), (int)C.dim1);
			
			
			
		#endif

	}
	C.isZeros = false;
	C.isIdentity = false;
#ifdef __COMPILE_TIMERS__
A.t_gemm.stop();
#endif
}



template<typename T>
void add(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C)
{   
#ifdef __COMPILE_TIMERS__
A.t_add_ext.start();
#endif
	if ((A.getNbCol() != B.getNbCol()) || (A.getNbRow() != B.getNbRow()) || (A.getNbRow() != C.getNbRow()) || (A.getNbCol() != C.getNbCol()))
	{
		handleError("LinAlgebra"," add : matrix dimension not equal");
		
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
template<typename T>	
T power_iteration(const  faust_mat<T> & A, const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag)
{	
	#ifdef __COMPILE_TIMERS__
		A.t_power_iteration.start();
	#endif


   const int nb_col = A.getNbCol();
   int i = 0;
   flag = 0;
	 
   if (nbr_iter_max <= 0)
      handleError("LinAlgebra "," power_iteration :  nbr_iter_max <= 0");
   if (nb_col != A.getNbRow())
      handleError("LinAlgebra "," power_iteration : faust_core<T> 1 must be a squared matrix"); 	
	 
   faust_vec<T> xk(nb_col);
   xk.setOnes();
   faust_vec<T> xk_norm(nb_col);
   T lambda_old=1.0;
   T lambda = 0.0;
   T alpha = 1.0;
   T beta = 0.0;
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
   return lambda;



 	
	 /*faust_unsigned_int nb_col = A.getNbCol();
	 faust_unsigned_int nb_row = A.getNbRow();
	 faust_unsigned_int i = 0;
	 faust_unsigned_int k;
	 bool do_continue=true;
	 T abs_eigen_value;
	 
	 bool stop_crit;
	 flag = 0;

	 if (nbr_iter_max <= 0)
		handleError("LinAlgebra","power_iteration : nbr_iter_max <= 0");
	 if (nb_col != nb_row)
		handleError("LinAlgebra","power_iteration : A must be square");
	 
	 faust_vec<T> xk(nb_col);
	 faust_vec<T> xk_pp(nb_col);
	 xk.setOnes();
	 xk.normalize();
	 
	 while(do_continue)
	 {	
		i++;
		gemv<T>(A,xk,xk_pp,1.,0,'N');
		abs_eigen_value = xk_pp.norm();
		//std::cout<<"current_norm : "<< abs_eigen_value<<std::endl;
		xk_pp.scalarMultiply(1/abs_eigen_value);
		
		stop_crit =true;
		k=0;
		while ((stop_crit) && (k<nb_col))
		{
			if (fabs(xk_pp(k) - xk(k))>threshold)
			{
				stop_crit = false;
			}
			k++;
		}
		
		if (stop_crit)
		{
			do_continue = false;
			flag = i;
			//std::cout<<"convergence : "<< i <<std::endl;

		}
		
		if (i >= nbr_iter_max)
		{	//std::cout<<"divergence"<<std::endl;
			do_continue = false;
			flag = -1;
		}
		
		
		xk = xk_pp;
      //std::cout << "i = " << i << " ; lambda=" << abs_eigen_value << std::endl;
	 }
	 //std::cout<<" flag :"<<flag<<std::endl;
	 #ifdef __COMPILE_TIMERS__
		A.t_power_iteration.stop();
	#endif
	//std::cout<<"flag inside power_it : "<<flag<<std::endl;
	//std::cout<<"threshold inside power_it : "<<threshold<<std::endl;
	//std::cout<<"max_it inside power_it : "<<nbr_iter_max<<std::endl;	
	 return abs_eigen_value;*/
	 
}




// non-member operators definitions
template<typename T>
faust_vec<T> operator*(const faust_mat<T>& M, const faust_vec<T>& v)
{
	faust_vec<T> vec(v);
	vec.multiplyLeft(M);
	return vec;
}

#endif
