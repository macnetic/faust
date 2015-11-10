#include <iostream>
#include "LinAlgebra.h"
#include "faust_mat.h"
#include "faust_vec.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "faust_spmat.h"
#include "faust_core.h"
#include "faust_exception.h"

	//////////FONCTION faust_mat<T> - faust_mat<T> ////////////////////

#ifdef __COMPILE_TIMERS__
	#include "faust_timer.h"
#endif

#ifdef __GEMM_WITH_OPENBLAS__
	#include "cblas.h"
#endif

// const char * core_algebra_name="faust_core<T>_algebra : ";

template<typename T>
T power_iteration(const  faust_core<T> & A, const int nbr_iter_max,T threshold, int & flag)
{	
	
	 int nb_col = A.getNbCol();
	 int nb_row = A.getNbRow();
	 int i = 0;
	 int k;
	 bool do_continue=true;
	 T abs_eigen_value;
	 
	 bool stop_crit;
	 flag = 0;

	 
	 if (nbr_iter_max <= 0)
	 {
		handleError("faust_core algebra "," power_iteration :  nbr_iter_max <= 0");
	 }
	 if (nb_col != nb_row)
	 {
		handleError("faust_core algebra "," power_iteration : faust_core<T> 1 must be a squared matrix"); 	
	 }
	 
	 
	 faust_vec<T> xk(nb_col);
	 faust_vec<T> xk_pp(nb_col);
	 xk.setOnes();
	 xk.normalize();
	 
	 while(do_continue)
	 {	
		i++;
		//gemv(A,xk,xk_pp,1.,0,'N');
		xk_pp = A*xk;
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
	 }
	 //std::cout<<" flag :"<<flag<<std::endl;

	/*std::cout<<"flag inside power_it : "<<flag<<std::endl;
	std::cout<<"threshold inside power_it : "<<threshold<<std::endl;
	std::cout<<"max_it inside power_it : "<<nbr_iter_max<<std::endl;*/		
	 return abs_eigen_value;
	 
}


// template<typename T>
// void multiply(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult)
 // {
	// int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB, nb_fact;

	// if  ((&(C.mat)) == (&(B.mat))) 
	// {
		// handleError("faust_core algebra "," multiply : C is the same object as B"); 		
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
				// handleError("faust_core algebra "," multiply :  dimension of faust_core<T> 1 and faust_spmat mismatch");	
			// }
		// }else
		// {
		
		
			// if (nbColOpB != nbRowOpA)
			// {

				// handleError("faust_core algebra "," multiply : dimension of faust_core<T> A and faust_spmat B mismatch");		
			// }
		// }
	// }else
	// {
		// handleWarning(" faust_core<T>_algebra : multiply : empty faust_core<T>");
	// }
	
	// if the faust_core<T> A is empty, it's considere as the identity, so C = equal B, it is useful into the algorithm palm4MSA, where the faust_core<T>s L and R can be empty	
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




template<typename T>
faust_vec<T> operator*(const faust_core<T>& f, const faust_vec<T>& v)
{
	faust_vec<T> vec(v);
	if (f.size() == 0)
		handleWarning("faust_core<T> algebra : operator* : empty faust_core<T>");
	
	for (int i=f.size()-1 ; i >= 0 ; i--)
		vec.multiplyLeft(f.data[i]);
	return vec;
}

template<typename T>
faust_mat<T> operator*(const faust_core<T>& f, const faust_mat<T>& M)
{
	faust_mat<T> A(M);
	if (f.size() == 0)
		handleWarning("faust_core<T> algebra : operator * : empty faust_core<T>");
	
	for (int i=f.size()-1 ; i >= 0 ; i--)
		A.multiplyLeft(f.data[i]);
	return A;
}



template<typename T>
void faust_solve(const faust_core<T> & A,faust_vec<T> & x, const faust_vec<T> & y)
{
	if (A.size() == 0)
	{
		handleError("faust_core algebra "," solve :  empty faust_core<T>");		
	}
	if (A.getNbRow() != y.size())
	{
		handleError("faust_core algebra "," solve :  dimension of the faust_core<T> A and size of y mismatch.");
	}
	faust_vec<T> ybis = y;
	
	for (int i=A.size()-1;i>=0;i--)
	{
		sp_solve(A.get_fact(i),x,ybis);
		ybis = x;
	}
	
	
	
}
