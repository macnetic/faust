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
 template<typename FPP>
 void Faust::multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult)
  {
	 int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB, nbFact;

	 if  ((&(C.mat)) == (&(B.mat)))
	 {
		 handleError("Faust::Transform algebra "," Faust::multiply : C is the same object as B");
	 }

	 nbFact = A.size();
	 if (nbFact != 0)
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


	 if (nbFact != 0)
	 {
		 if (typeMult == 'R')
		 {
			 if (nbColOpA != nbRowOpB)
			 {
				 handleError("Faust::Transform algebra "," Faust::multiply :  dimension of Faust::Transform<FPP,Cpu> 1 and Faust::MatSparse mismatch");
			 }
		 }else
		 {


			 if (nbColOpB != nbRowOpA)
			 {

				 handleError("Faust::Transform algebra "," Faust::multiply : dimension of Faust::Transform<FPP,Cpu> A and Faust::MatSparse B mismatch");
			 }
		 }
	 }else
	 {
		 handleWarning(" Faust::Transform<FPP,Cpu>_algebra : Faust::multiply : empty Faust::Transform<FPP,Cpu>");
	 }

	 // if the Faust::Transform<FPP,Cpu> A is empty, it's considere as the identity, so C = equal B, it is useful into the algorithm Faust::Palm4MSA, where the Faust::Transform<FPP,Cpu>s L and R can be empty
	 C = B;
	 C.scalarMultiply(alpha);
	 C.resize(nbRowOpB,nbColOpB);
	 if (nbFact != 0)
	 {
		 if (typeA == 'T')
		 {
			 if(typeMult == 'R')
			 {
				 for (int i=0 ; i<nbFact ; i++)
				 {
					 C.mat = A.data[i].mat.transpose() * C.mat;
				 }
				 C.resize(nbRowOpA,nbColOpB);

			 }else
			 {
				 for (int i=nbFact-1 ; i>=0 ; i--)
				 {
				 C.mat = C.mat * A.data[i].mat.transpose();
				 }
				 C.resize(nbRowOpB,nbColOpA);
			 }


		 }else
		 {
			 if(typeMult == 'R')
			 {
				 for (int i=nbFact-1 ; i>=0 ; i--)
				 {
					 C.mat = A.data[i].mat * C.mat;
				 }
				 C.resize(nbRowOpA,nbColOpB);
			 }else
			 {
				 for (int i=0 ; i<nbFact ; i++)
				 {
					 C.mat = C.mat*A.data[i].mat;
				 }
				 C.resize(nbRowOpB,nbColOpA);
			 }
		 }
	 }


  }
//////////// modif AL AL



template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::Vect<FPP,Cpu>& v)
{
	Faust::Vect<FPP,Cpu> vec(v);
	if (f.size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> algebra : operator* : empty Faust::Transform<FPP,Cpu>");

	for (int i=f.size()-1 ; i >= 0 ; i--)
		vec.multiplyLeft(f.data[i]);
	return vec;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::operator*(const Faust::Transform<FPP,Cpu>& f, const Faust::MatDense<FPP,Cpu>& M)
{
	Faust::MatDense<FPP,Cpu> A(M);
	if (f.size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> algebra : operator * : empty Faust::Transform<FPP,Cpu>");

	for (int i=f.size()-1 ; i >= 0 ; i--)
		A.multiplyLeft(f.data[i]);
	return A;
}

#endif
