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

	//////////FONCTION faust_mat - faust_mat ////////////////////

#ifdef __COMPILE_TIMERS__
	#include "faust_timer.h"
#endif

#ifdef __GEMM_WITH_OPENBLAS__
	#include "cblas.h"
#endif

faust_real power_iteration(const  faust_core & A, const int nbr_iter_max,faust_real threshold, int & flag)
{	
	
	 int nb_col = A.getNbCol();
	 int nb_row = A.getNbRow();
	 int i = 0;
	 int k;
	 bool do_continue=true;
	 faust_real abs_eigen_value;
	 
	 bool stop_crit;
	 flag = 0;

	 
	 if (nbr_iter_max <= 0)
	 {
		handleError(" Faust_core_algebra : power_iteration :  nbr_iter_max <= 0");
	 }
	 if (nb_col != nb_row)
	 {
		handleError(" Faust_core_algebra : power_iteration : faust_core 1 must be a squared matrix"); 	
	 }
	 
	 
	 faust_vec xk(nb_col);
	 faust_vec xk_pp(nb_col);
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

 void multiply(const faust_core & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, char typeA, char typeMult)
 {
	int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB, nb_fact;

	if  ((&(C.mat)) == (&(B.mat))) 
	{
		handleError(" Faust_core_algebra : multiply : C is the same object as B"); 		
	}
	
	nb_fact = A.size();
	if (nb_fact != 0)
	{	
		if (typeA == 'T')
		{
			nbRowOpA = A.getNbCol();
			nbColOpA = A.getNbRow();
		}else
		{
			nbRowOpA = A.getNbRow();
			nbColOpA = A.getNbCol();
		}
	}

		nbRowOpB = B.getNbRow();
		nbColOpB = B.getNbCol();	
	
	
	if (nb_fact != 0)
	{
		if (typeMult == 'R')
		{	
			if (nbColOpA != nbRowOpB)
			{
				handleError(" Faust_core_algebra : multiply :  dimension of faust_core 1 and faust_spmat mismatch");	
			}
		}else
		{
		
		
			if (nbColOpB != nbRowOpA)
			{

				handleError(" Faust_core_algebra : multiply : dimension of faust_core A and faust_spmat B mismatch");		
			}
		}
	}else
	{
		handleWarning(" Faust_core_algebra : multiply : empty faust_core");
	}
	
	// if the faust_core A is empty, it's considere as the identity, so C = equal B, it is useful into the algorithm palm4MSA, where the faust_cores L and R can be empty	
	C = B;
	C.scalarMultiply(alpha);
	C.resize(nbRowOpB,nbColOpB);
	if (nb_fact != 0)
	{
		if (typeA == 'T')
		{
			if(typeMult == 'R')
			{
				for (int i=0 ; i<nb_fact ; i++)
				{	
					C.mat = A.data[i].mat.transpose() * C.mat;
				}
				C.resize(nbRowOpA,nbColOpB);	
				
			}else
			{
				for (int i=nb_fact-1 ; i>=0 ; i--)
				{
				C.mat = C.mat * A.data[i].mat.transpose();
				}
				C.resize(nbRowOpB,nbColOpA);
			}
				
				
		}else
		{	
			if(typeMult == 'R')
			{
				for (int i=nb_fact-1 ; i>=0 ; i--)
				{	
					C.mat = A.data[i].mat * C.mat;
				}
				C.resize(nbRowOpA,nbColOpB);
			}else
			{
				for (int i=0 ; i<nb_fact ; i++)
				{	
					C.mat = C.mat*A.data[i].mat;
				}
				C.resize(nbRowOpB,nbColOpA);
			}
		}	
	}
	
	
 }





faust_vec operator*(const faust_core& f, const faust_vec& v)
{
	faust_vec vec(v);
	if (f.size() == 0)
		handleWarning("faust_core algebra : operator* : empty faust_core");
	
	for (int i=f.size()-1 ; i >= 0 ; i--)
		vec.multiplyLeft(f.data[i]);
	return vec;
}

faust_mat operator*(const faust_core& f, const faust_mat& M)
{
	faust_mat A(M);
	if (f.size() == 0)
		handleWarning("faust_core algebra : operator * : empty faust_core");
	
	for (int i=f.size()-1 ; i >= 0 ; i--)
		A.multiplyLeft(f.data[i]);
	return A;
}
