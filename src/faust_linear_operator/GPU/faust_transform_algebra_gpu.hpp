#ifndef __FAUST_TOOLS_ALGEBRA_GPU_HPP__
#define __FAUST_TOOLS_ALGEBRA_GPU_HPP__

#include <iostream>
#include "faust_linear_algebra_gpu.h"
#include "faust_MatDense_gpu.h"
#include "faust_Vect_gpu.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "faust_MatSparse_gpu.h"
#include "faust_Transform_gpu.h"
#include "faust_exception.h"

	//////////FONCTION faust_cu_mat<T> - faust_cu_mat<T> ////////////////////

#ifdef __COMPILE_TIMERS__
	#include "faust_Timer.h"
#endif

#include "cublas_v2.h"

// const char * core_algebra_name="faust_core_cu<T>_algebra : ";
/*
template<typename T>
T power_iteration(const  faust_core_cu<T> & cu_A, const int nbr_iter_max, const T threshold, int & flag, cublasHandle_t cublasHandle)
{
   const int nb_col = cu_A.getNbCol();
   int i = 0;
   flag = 0;

   if (nbr_iter_max <= 0)
      handleError("faust_core_cu algebra "," power_iteration :  nbr_iter_max <= 0");
   if (nb_col != cu_A.getNbRow())
      handleError("faust_core_cu algebra "," power_iteration : faust_core_cu<T> 1 must be a squared matrix");

   faust_cu_vec<T> cu_xk(nb_col);
   cu_xk.setOnes();
   faust_cu_vec<T> cu_xk_norm(nb_col);
   T lambda_old=1.0;
   T lambda = 0.0;
   while(fabs(lambda_old-lambda)>threshold && i<nbr_iter_max)
   {
      i++;
      lambda_old = lambda;
      cu_xk_norm = cu_xk;
      cu_xk_norm.normalize();
      gemv(cu_A, cu_xk_norm, cu_xk, cublasHandle);
      lambda = cu_xk_norm.dot(cu_xk);
   }
   flag = (i<nbr_iter_max)?i:-1;
   return lambda;
}
*/

// template<typename T>
// void multiply(const faust_core_cu<T> & A, const faust_cu_mat<T> & B, faust_cu_mat<T> & C,const T & alpha, char typeA, char typeMult)
 // {
	// int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB, nb_fact;

	// if  ((&(C.mat)) == (&(B.mat)))
	// {
		// handleError("faust_core_cu algebra "," multiply : C is the same object as B");
	// }

	// nb_fact = A.size();
	// if (nb_fact != 0)
	// {
		// if (typeA == 'T')
		// {
			// nbRowOpA = A.getNbCol();
			// nbColOpA = A.getNbRow();
		// }else
		// {
			// nbRowOpA = A.getNbRow();
			// nbColOpA = A.getNbCol();
		// }
	// }

		// nbRowOpB = B.getNbRow();
		// nbColOpB = B.getNbCol();


	// if (nb_fact != 0)
	// {
		// if (typeMult == 'R')
		// {
			// if (nbColOpA != nbRowOpB)
			// {
				// handleError("faust_core_cu algebra "," multiply :  dimension of faust_core_cu<T> 1 and faust_spmat mismatch");
			// }
		// }else
		// {


			// if (nbColOpB != nbRowOpA)
			// {

				// handleError("faust_core_cu algebra "," multiply : dimension of faust_core_cu<T> A and faust_spmat B mismatch");
			// }
		// }
	// }else
	// {
		// handleWarning(" faust_core_cu<T>_algebra : multiply : empty faust_core_cu<T>");
	// }

	// if the faust_core_cu<T> A is empty, it's considere as the identity, so C = equal B, it is useful into the algorithm Palm4MSA, where the faust_core_cu<T>s L and R can be empty
	// C = B;
	// C.scalarMultiply(alpha);
	// C.resize(nbRowOpB,nbColOpB);
	// if (nb_fact != 0)
	// {
		// if (typeA == 'T')
		// {
			// if(typeMult == 'R')
			// {
				// for (int i=0 ; i<nb_fact ; i++)
				// {
					// C.mat = A.data[i].mat.transpose() * C.mat;
				// }
				// C.resize(nbRowOpA,nbColOpB);

			// }else
			// {
				// for (int i=nb_fact-1 ; i>=0 ; i--)
				// {
				// C.mat = C.mat * A.data[i].mat.transpose();
				// }
				// C.resize(nbRowOpB,nbColOpA);
			// }


		// }else
		// {
			// if(typeMult == 'R')
			// {
				// for (int i=nb_fact-1 ; i>=0 ; i--)
				// {
					// C.mat = A.data[i].mat * C.mat;
				// }
				// C.resize(nbRowOpA,nbColOpB);
			// }else
			// {
				// for (int i=0 ; i<nb_fact ; i++)
				// {
					// C.mat = C.mat*A.data[i].mat;
				// }
				// C.resize(nbRowOpB,nbColOpA);
			// }
		// }
	// }


 // }




//template<typename T>
//faust_cu_vec<T> operator*(const faust_core_cu<T>& f, const faust_cu_vec<T>& v)
//{
//	faust_cu_vec<T> vec(v);
//	if (f.size() == 0)
//		handleWarning("faust_core_cu<T> algebra : operator* : empty faust_core_cu<T>");
//
//	for (int i=f.size()-1 ; i >= 0 ; i--)
//		vec.multiplyLeft(f.data[i]);
//	return vec;
//}

//template<typename T>
//faust_cu_mat<T> operator*(const faust_core_cu<T>& f, const faust_cu_mat<T>& M)
//{
//	faust_cu_mat<T> A(M);
//	if (f.size() == 0)
//		handleWarning("faust_core_cu<T> algebra : operator * : empty faust_core_cu<T>");

//	for (int i=f.size()-1 ; i >= 0 ; i--)
//		A.multiplyLeft(f.data[i]);
//	return A;
//}


#endif
