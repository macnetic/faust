#ifndef __FAUST_Transform_ALGEBRA_HPP__
#define __FAUST_Transform_ALGEBRA_HPP__

#include <iostream>
#include "faust_linear_algebra.h"
#include "faust_MatDense.h"
#include "faust_Vect.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "faust_MatSparse.h"
#include "faust_Transform.h"
#include "faust_exception.h"

	//////////FONCTION Faust::MatDense<FPP,Cpu> - Faust::MatDense<FPP,Cpu> ////////////////////

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif



// const char * core_algebra_name="Faust::Transform<FPP,Cpu>_algebra : ";
template<typename FPP>
FPP Faust::power_iteration(const  Faust::Transform<FPP,Cpu> & A, const int nbr_iter_max,FPP threshold, int & flag)
{


 const int nb_col = A.getNbCol();
   int i = 0;
   flag = 0;

   if (nbr_iter_max <= 0)
      handleError("faust_transform_algebra "," power_iteration :  nbr_iter_max <= 0");
   if (nb_col != A.getNbRow())
      handleError("faust_transform_algebra "," power_iteration : Faust::Transform_gpu<FPP,Cpu> 1 must be a squared matrix");

   Faust::Vect<FPP,Cpu> xk(nb_col);
   xk.setOnes();
   Faust::Vect<FPP,Cpu> xk_norm(nb_col);
   FPP lambda_old=1.0;
   FPP lambda = 0.0;
   while(fabs(lambda_old-lambda)>threshold && i<nbr_iter_max)
   {
      i++;
      lambda_old = lambda;
      xk_norm = xk;
      xk_norm.normalize();
      xk = A*xk_norm;
      lambda = xk_norm.dot(xk);
   }
   flag = (i<nbr_iter_max)?i:-1;
   return lambda;










}

//////////// modif AL AL
 
//////////// modif AL AL

template<typename FPP>
#ifdef __COMPILE_TIMERS__
	Faust::Vect<FPP,Cpu> Faust::operator*(Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v)
#else		
	Faust::Vect<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu> & v)
#endif
{

	return f.multiply(v);
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M)
{
	return f.multiply(M);
}

#endif
