#include "palm4MSA.h"

#include <iostream>
#include "faust_params.h"
#include "faust_params_palm.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_real.h"
#include "faust_constraint_int.h"
#include "faust_mat.h"
#include "LinAlgebra.h"
#include "prox.h"

#ifdef __COMPILE_TIMERS__
#include "faust_timer.h"
#endif


//#include"faust_init_from_matio_mat.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#define __SP setprecision(20)<<

using namespace std;

palm4MSA::palm4MSA(const faust_params& params_, const bool isGlobal_) :
   data(params_.data),
   lambda(params_.init_lambda),
   nb_fact(0),
   S(params_.init_fact),
   RorL(vector<faust_mat>(2)),
   ind_fact(0),
   ind_ite(-1),
   lipschitz_multiplicator(1.001),
   verbose(params_.isVerbose),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isLambdaComputed(params_.isLambdaComputed),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false),
   isGlobal(isGlobal_),
   isInit(false)
{
   if(isGlobal)
      stop_crit = stopping_criterion(params_.stop_crit_global);
   else
      stop_crit = stopping_criterion(params_.stop_crit_2facts);
}

palm4MSA::palm4MSA(const faust_params_palm& params_palm_, const bool isGlobal_) :
   stop_crit(params_palm_.stop_crit),
   data(params_palm_.data),
   lambda(params_palm_.init_lambda),
   nb_fact(params_palm_.nb_fact),
   S(params_palm_.init_fact),
   RorL(vector<faust_mat>(2)),
   LorR(faust_mat(params_palm_.init_fact[0].getNbRow())),
   const_vec(params_palm_.cons),
   ind_fact(0),
   ind_ite(-1),
   lipschitz_multiplicator(1.001),
   verbose(params_palm_.isVerbose),
   isUpdateWayR2L(params_palm_.isUpdateWayR2L),
   isLambdaComputed(params_palm_.isLambdaComputed),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false),
   isGlobal(isGlobal_)
{
   check_constraint_validity();

}

void palm4MSA::compute_facts()
{
	while (do_continue())
	{	
		next_step();
		//cout<<"lambda : "<<lambda<<endl;
	}
	//cout<<"lambda_computed : "<<isLambdaComputed<<endl;
	if (!isLambdaComputed)
	{	
		isLastFact=true;
		//cout<<"lambda before : "<<lambda<<endl;
		compute_last_update();
		//cout<<"final lambda : "<<lambda<<endl;
		
	}
}

void palm4MSA::compute_projection()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.start();
t_local_compute_projection.start();
#endif

   if (const_vec[ind_fact]->getConstraintType() == CONSTRAINT_NAME_CONST)
   {	
		#ifdef __COMPILE_TIMERS__
			nb_call_prox_const++;
			t_prox_const.start();
			
		#endif
         const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
         S[ind_fact] = const_mat->getParameter();
		 #ifdef __COMPILE_TIMERS__
			t_prox_const.stop();
		#endif
   }
   else
   {
      //faust_mat matrix2project(S[ind_fact]);
      S[ind_fact] -= grad_over_c;
      switch (const_vec[ind_fact]->getConstraintType())
      {
         case CONSTRAINT_NAME_SP:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_sp++;
				t_prox_sp.start();
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
			if (isLambdaComputed)
			{
				prox_sp(S[ind_fact], const_int->getParameter());
			}
			else
			{	
				prox_sp_normfree(S[ind_fact], const_int->getParameter());
			}
			#ifdef __COMPILE_TIMERS__
			t_prox_sp.stop();
			#endif
			
			
         }
         break;

         case CONSTRAINT_NAME_SPCOL:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_spcol++;
				t_prox_spcol.start();
				
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
			if (isLambdaComputed)
			{
            prox_spcol(S[ind_fact], const_int->getParameter());
			}else
			{
				prox_spcol_normfree(S[ind_fact], const_int->getParameter());	
			}
			#ifdef __COMPILE_TIMERS__
				t_prox_spcol.stop();
			#endif
         }
         break;

         case CONSTRAINT_NAME_SPLIN:
         {	

			

			#ifdef __COMPILE_TIMERS__
			nb_call_prox_splin++;
			t_prox_splin.start();
			#endif
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);

			if (isLambdaComputed)
			{
				prox_splin(S[ind_fact], const_int->getParameter());
			}else
			{
				prox_splin_normfree(S[ind_fact], const_int->getParameter());
			}
			#ifdef __COMPILE_TIMERS__
				t_prox_splin.stop();
			#endif
			
         }
         break;

         case CONSTRAINT_NAME_NORMCOL:
         {	
			#ifdef __COMPILE_TIMERS__
				nb_call_prox_normcol++;
				t_prox_normcol.start();
			#endif
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            prox_normcol(S[ind_fact], const_real->getParameter());
			#ifdef __COMPILE_TIMERS__
				t_prox_normcol.stop();
			#endif
         }
         break;

         case CONSTRAINT_NAME_SPLINCOL:
         {	

            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            //prox_splincol(S[ind_fact], const_real->getParameter());

         }
         break;

         case CONSTRAINT_NAME_L0PEN:
         {	

            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l0pen(S[ind_fact], sqrt(2*const_int->getParameter()/c));

         }
         break;

         case CONSTRAINT_NAME_L1PEN:
         {	

            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l1pen(S[ind_fact], const_int->getParameter()/c);

         }
         break;

         case CONSTRAINT_NAME_WAV:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_wav(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SP_POS:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
			if (isLambdaComputed)
			{
				prox_sp_pos(S[ind_fact], const_int->getParameter());
			}
			else
			{	
				prox_sp_pos_normfree(S[ind_fact], const_int->getParameter());
			}
         }
         break;

         case CONSTRAINT_NAME_BLKDIAG:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLIN_TEST:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SUPP:
         {
            const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_mat->getParameter());
         }
         break;

         case CONSTRAINT_NAME_NORMLIN:
         {
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            prox_normlin(S[ind_fact], const_real->getParameter());
         }
         break;

         case CONSTRAINT_NAME_TOEPLITZ:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         default:
            cerr << "error in palm4MSA::compute_projection : unknown name of constraint" << endl;
            exit(EXIT_FAILURE);

      }
   }
   isProjectionComputed = true;
#ifdef __COMPILE_TIMERS__
t_global_compute_projection.stop();
t_local_compute_projection.stop();
#endif
}


void palm4MSA::compute_grad_over_c()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_grad_over_c.start();
t_local_compute_grad_over_c.start();
#endif
   if(!isCComputed) 
   {
      cerr << "c must be set before computing grad/c" << endl;
      exit(EXIT_FAILURE);
   }

// There are 4 ways to compute gradient :
// (0) : lambda*(L'*(lambda*(L*S)*R - X))*R' : complexity = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2;
// (1) : lambda*L'*((lambda*(L*S)*R - X)*R') : complexity = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2;
// (2) : lambda*(L'*(lambda*L*(S*R) - X))*R' : complexity = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2;
// (3) : lambda*L'*((lambda*L*(S*R) - X)*R') : complexity = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2;
//  with L of size L1xL2
//       S of size L2xS2
//       R of size S2xR2

   unsigned long long int L1, L2, R2, S2;
   if (!isUpdateWayR2L)
   {
      L1 = (unsigned long long int) LorR.getNbRow();
      L2 = (unsigned long long int) LorR.getNbCol();
      R2 = (unsigned long long int) RorL[ind_fact].getNbCol();
   }
   else
   {
      L1 = (unsigned long long int) RorL[ind_fact].getNbRow();
      L2 = (unsigned long long int) RorL[ind_fact].getNbCol(); 
      R2 = (unsigned long long int) LorR.getNbCol();
   }
   S2 = (unsigned long long int) S[ind_fact].getNbCol();
   vector<unsigned long long int > complexity(4,0);
   complexity[0] = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2;
   complexity[1] = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2;
   complexity[2] = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2;
   complexity[3] = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2;

   int idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));


   error = data;
   faust_mat tmp1, tmp2, tmp3, tmp4;

   if (idx==0 || idx==1) // computing L*S first, then (L*S)*R
   {
      if (!isUpdateWayR2L)
      {
         // tmp1 = L*S
         multiply(LorR, S[ind_fact], tmp1);
         // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
         gemm(tmp1, RorL[ind_fact], error, lambda, -1.0, 'N', 'N');
      }
      else
      {
         // tmp1 = L*S
         multiply(RorL[ind_fact], S[ind_fact], tmp1);
         // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
         gemm(tmp1, LorR, error, lambda, -1.0, 'N', 'N');
      }
   }
   else // computing S*R first, then L*(S*R)
   {
      if (!isUpdateWayR2L)
      {
         // tmp1 = S*R
         multiply(S[ind_fact], RorL[ind_fact], tmp1);
         // error = lambda*L*tmp1 - error (= lambda*L*S*R - data )
         gemm(LorR, tmp1, error, lambda, -1.0, 'N', 'N');
      }
      else
      {
         // tmp1 = S*R
         multiply(S[ind_fact], LorR, tmp1);
         // error = lambda*L*tmp1 - error (= lambda*L*S*R - data )
         gemm(RorL[ind_fact], tmp1, error, lambda, -1.0, 'N', 'N');
      }
   }



   if (idx==0 || idx==2) // computing L'*error first, then (L'*error)*R'
   {
      if (!isUpdateWayR2L)
      {
         // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
         gemm(LorR, error, tmp3, lambda, 0.0, 'T', 'N');
         // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
         gemm(tmp3, RorL[ind_fact], grad_over_c, 1.0/c, 0.0,'N','T');
      }
      else
      {
         // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
         gemm(RorL[ind_fact], error, tmp3, lambda, 0.0, 'T', 'N');
         // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
         gemm(tmp3, LorR, grad_over_c, 1.0/c, 0.0,'N','T');
      }
   }
   else // computing error*R' first, then L'*(error*R')
   {
      if (!isUpdateWayR2L)
      {
         // tmp3 = lambda*error*R' (= lambda*(lambda*L*S*R - data) * R' )
         gemm(error, RorL[ind_fact], tmp3, lambda, 0.0, 'N', 'T');
         // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * lambda*(lambda*L*S*R - data) * R' )
         gemm(LorR, tmp3, grad_over_c,1.0/c, 0.0,'T','N');
      }
      else
      {
         // tmp3 = lambda*error*R' (= lambda * (lambda*L*S*R - data) * R' )
         gemm(error, LorR, tmp3, lambda, 0.0, 'N', 'T');
         // grad_over_c = 1/c*L'*tmp3 (= 1/c*L' * lambda*(lambda*L*S*R - data) * R' )
         gemm(RorL[ind_fact], tmp3, grad_over_c, 1.0/c, 0.0,'T','N');
      }

   }


   isGradComputed = true;

#ifdef __COMPILE_TIMERS__
t_global_compute_grad_over_c.stop();
t_local_compute_grad_over_c.stop();
#endif
}
void palm4MSA::compute_last_update()
{
	/*if (do_continue())
	{
		cerr << "error in palm4MSA::compute_last_update : computation of lambda must be done at the end of the main iteration " << endl;
		exit(EXIT_FAILURE);
	}*/
	
	if (isLambdaComputed)
	{
		cerr << "error in palm4MSA::compute_last_update : computation of lambda at the end whereas isLambdaComputed means it must be computed at each step" << endl;
		exit(EXIT_FAILURE);
	}
	
	faust_real scaling = 1;
	for (int i=0;i<nb_fact;i++)
	{
		
		faust_real current_norm;
		switch (const_vec[i]->getConstraintType())
		{
			case CONSTRAINT_NAME_NORMCOL:
			{}
			break;
			case CONSTRAINT_NAME_NORMLIN: 
			{}
			default:
			{
				current_norm=1/(S[i].norm());
				scaling = scaling*current_norm;
				S[i].scalarMultiply(current_norm);
			}
			
		}
	}
	LorR.scalarMultiply(scaling);
	//LorR = S[0];
	//for (int i=1;i<nb_fact;i++)LorR*=S[i];		
	//cout<<" scaling : "<<scaling<<endl; 
	compute_lambda();
	
}

void palm4MSA::compute_lambda()
{
#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.start();
t_local_compute_lambda.start();
#endif

   if (!isLastFact)
   {
      cerr << "error in palm4MSA::compute_lambda : computation of lambda must be done at the end of the iteration through the number of factors" << endl;
      exit(EXIT_FAILURE);
   }

   // As LorR has also been updated at the end of the last iteration over the facts, LorR matches X_hat, which the product of all factors, including the last one.
   // Xt_Xhat = data'*X_hat
   faust_mat Xt_Xhat;
   gemm(data, LorR, Xt_Xhat, 1.0, 0.0, 'T','N');

   // Xhatt_Xhat = X_hat'*X_hat
   faust_mat Xhatt_Xhat;
   gemm(LorR, LorR, Xhatt_Xhat, 1.0, 0.0, 'T','N');


   lambda = Xt_Xhat.trace()/Xhatt_Xhat.trace();

   //cout<<"lambda : "<<lambda<<endl;
   //cout<<__SP lambda<<endl;

#ifdef __COMPILE_TIMERS__
t_global_compute_lambda.stop();
t_local_compute_lambda.stop();
#endif
}


void palm4MSA::update_R()
{
#ifdef __COMPILE_TIMERS__
t_global_update_R.start();
t_local_update_R.start();
#endif

   // R[nb_fact-1] est initialise a l'identite lors de la creation de l'objet palm4MSA et n'est pas cense changer
   if (!isUpdateWayR2L)
   {
      RorL.resize(nb_fact);
      RorL[nb_fact-1].resize(const_vec[nb_fact-1]->getCols());
      RorL[nb_fact-1].setEyes();

      for (int i=nb_fact-2 ; i>-1 ; i--)
         //  R[i] = S[i+1] * R[i+1]
         multiply(S[i+1], RorL[i+1], RorL[i]);
   }
   else
   {
      if(!isProjectionComputed)
      {
         cerr << "Projection must be computed before updating L" << endl;
         exit(EXIT_FAILURE);
      }
      LorR.multiplyLeft(S[ind_fact]);
   }
#ifdef __COMPILE_TIMERS__
t_global_update_R.stop();
t_local_update_R.stop();
#endif
}


void palm4MSA::update_L()
{
#ifdef __COMPILE_TIMERS__
t_global_update_L.start();
t_local_update_L.start();
#endif

   if(!isUpdateWayR2L)
   {
      if(!isProjectionComputed)
      {
         cerr << "Projection must be computed before updating L" << endl;
         exit(EXIT_FAILURE);
      }
      LorR *= S[ind_fact];
   }
   else
   {
      RorL.resize(nb_fact);
      RorL[0].resize(const_vec[0]->getRows());
      RorL[0].setEyes();
      for (int i=0 ; i<nb_fact-1 ; i++)
         //  R[i] = S[i+1] * R[i+1]
         multiply(RorL[i] , S[i], RorL[i+1]);
   }
#ifdef __COMPILE_TIMERS__
t_global_update_L.stop();
t_local_update_L.stop();
#endif
}

void palm4MSA::check_constraint_validity()
{
#ifdef __COMPILE_TIMERS__
t_global_check.start();
t_local_check.start();
#endif

   if (nb_fact != S.size())
   {
      cerr << "Error in palm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << endl;
      exit(EXIT_FAILURE);
   }
#ifdef __COMPILE_TIMERS__
t_global_check.stop();
t_local_check.stop();
#endif
}


#ifdef __PAS_FIXE__
// palm4MSA::compute_c() has been defined as an inline method in palm4MSA.h
#else
void palm4MSA::compute_c()
{
#ifdef __COMPILE_TIMERS__
	t_global_compute_c.start();
#endif

   //std::cout<<"calcul pas : "<<std::endl;
   // faust_real nL=LorR.spectralNorm();
   //faust_real nR=RorL[ind_fact].spectralNorm();
   //c=lipschitz_multiplicator*nR*nR*nL*nL*lambda*lambda;


   int flag1,flag2;
   
   int nbr_iter = 10000;
   faust_real threshold = 1e-16;
   faust_real nL1=LorR.spectralNorm(nbr_iter,threshold,flag1);
   faust_real nR1=RorL[ind_fact].spectralNorm(nbr_iter,threshold,flag2);
    c=lipschitz_multiplicator*nR1*nR1*nL1*nL1*lambda*lambda;

   //std::cout<<" nL : "<<nL <<" nL1 : "<<nL1<<" flag : "<<flag1<<std::endl;
   //std::cout<<" nR : "<<nR <<" nR1 : "<<nR1<<" flag : "<<flag2<<std::endl;
   //std::cout<<" c : "<< c <<" c1 : "<<c1<<std::endl<<std::endl;
   //std::cout<<" c : "<< c <<std::endl;


   isCComputed = true;
   
   #ifdef __COMPILE_TIMERS__
	t_global_compute_c.stop();
	#endif	
}
#endif


void palm4MSA::init_fact(int nb_facts_)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact.start();
t_local_init_fact.start();
#endif

  /*if(isInit && isGlobal)
  {
     cerr << "Error in palm4MSA::init_fact : global factorization has already been initialized" <<endl;
     exit(EXIT_FAILURE);
  }
  else if (isGlobal)
  {
     nb_fact = 1;
     S.resize(nb_fact);
     isInit = true;
     return;
  }*/



// if !isGlobal

  if(!isConstraintSet)
  {
     cerr << "Error in palm4MSA::init_fact : constrainst must be set before calling init_fact" << endl;
     exit(EXIT_FAILURE);
  }
   nb_fact = nb_facts_;
   S.resize(nb_fact);
   if (!isUpdateWayR2L)
   {
      S[0].resize(const_vec[0]->getRows(), const_vec[0]->getCols());
      S[0].setZeros();
      for (int i=1 ; i<nb_fact ; i++)
      {
         S[i].resize(const_vec[i]->getRows(), const_vec[i]->getCols());
         S[i].setEyes();   
      }   
   }
   else
   {
      for (int i=0 ; i<nb_fact-1 ; i++)
      {
         S[i].resize(const_vec[i]->getRows(), const_vec[i]->getCols());
         S[i].setEyes();    
      } 
      S[nb_fact-1].resize(const_vec[nb_fact-1]->getRows(), const_vec[nb_fact-1]->getCols());
      S[nb_fact-1].setZeros();   
   } 

    
#ifdef __COMPILE_TIMERS__
t_global_init_fact.stop();
t_local_init_fact.stop();
#endif
}

void palm4MSA::next_step()
{
#ifdef __COMPILE_TIMERS__
t_global_next_step.start();
t_local_next_step.start();
#endif


   check_constraint_validity();
   // resizing L or R 
   if(!isUpdateWayR2L)
   {
      LorR.resize(const_vec[0]->getRows());
      LorR.setEyes();
      update_R();
   }
   else
   {
      LorR.resize(const_vec[nb_fact-1]->getCols());
      LorR.setEyes();
      update_L();
   }

   int* ind_ptr = new int[nb_fact];
   for (int j=0 ; j<nb_fact ; j++)
      if (!isUpdateWayR2L)
         ind_ptr[j] = j;
      else
         ind_ptr[j] = nb_fact-1-j;
     
   for (int j=0 ; j<nb_fact ; j++)
   {
      if (j == nb_fact-1)
         isLastFact = true;
      else
         isLastFact = false;

      ind_fact = ind_ptr[j];

      isCComputed = false;
      isGradComputed = false;
      isProjectionComputed = false;
      
      compute_c();
      compute_grad_over_c();
      compute_projection();
  
 
      if(!isUpdateWayR2L)
         update_L();
      else
         update_R();

   }
	if (isLambdaComputed)
	{	
		compute_lambda();

	}
   if (verbose)
      cout << "Iter " << ind_ite << ", RMSE=" << get_RMSE() << endl;


   delete[] ind_ptr;
   ind_ptr = NULL;

//cout<<"lambda : "<< lambda<< endl;   
#ifdef __COMPILE_TIMERS__
t_global_next_step.stop();
t_local_next_step.stop();
#endif
}

void palm4MSA::init_fact_from_palm(const palm4MSA& palm2, bool isFactSideLeft)
{
#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.start();
t_local_init_fact_from_palm.start();
#endif


   if (palm2.nb_fact != 2)
   {
      cerr << "argument palm2 must contain 2 factors." << endl;
      exit(EXIT_FAILURE);
   }

  if(!isConstraintSet)
  {
     cerr << "Error in palm4MSA::init_fact_from_palm : constrainst must be set before calling init_fact_from_palm" << endl;
     exit(EXIT_FAILURE);
  }

   if(isFactSideLeft) 
   {
      S.insert(S.begin(), palm2.S[0]);
      S[1] = palm2.S[1];
   }
   else
   {
      S[S.size()-1] = palm2.S[0];
      S.push_back(palm2.S[1]);
   }
   nb_fact++;

   check_constraint_validity();

#ifdef __COMPILE_TIMERS__
t_global_init_fact_from_palm.stop();
t_local_init_fact_from_palm.stop();
#endif
}


#ifdef __COMPILE_TIMERS__

faust_timer palm4MSA::t_global_compute_projection;
faust_timer palm4MSA::t_global_compute_grad_over_c;
faust_timer palm4MSA::t_global_compute_c;
faust_timer palm4MSA::t_global_compute_lambda;
faust_timer palm4MSA::t_global_update_R;
faust_timer palm4MSA::t_global_update_L;
faust_timer palm4MSA::t_global_check;
faust_timer palm4MSA::t_global_init_fact;
faust_timer palm4MSA::t_global_next_step;
faust_timer palm4MSA::t_global_init_fact_from_palm;


faust_timer palm4MSA::t_prox_const;
faust_timer palm4MSA::t_prox_sp;
faust_timer palm4MSA::t_prox_spcol;
faust_timer palm4MSA::t_prox_splin;
faust_timer palm4MSA::t_prox_normcol;

int palm4MSA::nb_call_prox_const;
int palm4MSA::nb_call_prox_sp;
int palm4MSA::nb_call_prox_spcol;
int palm4MSA::nb_call_prox_splin;
int palm4MSA::nb_call_prox_normcol;





void palm4MSA::init_local_timers()
{
t_local_compute_projection.reset();
t_local_compute_grad_over_c.reset();
t_local_compute_lambda.reset();
t_local_update_R.reset();
t_local_update_L.reset();
t_local_check.reset();
t_local_init_fact.reset();
t_local_next_step.reset();
t_local_init_fact_from_palm.reset();
}

void palm4MSA::print_global_timers()const
{
   cout << "timers in palm4MSA : " << endl;
   cout << "t_global_next_step           = " << t_global_next_step.get_time()           << " s for "<< t_global_next_step.get_nb_call()           << " calls" << endl;
   cout << "t grad + updateL + updateR  = " << t_global_compute_grad_over_c.get_time() + t_global_update_L.get_time() + t_global_update_R.get_time()  << " s for "<< t_global_compute_grad_over_c.get_nb_call()            << " calls of grad" << endl;
     cout << "t_global_compute_c = " << t_global_compute_c.get_time() << " s for "<< t_global_compute_c.get_nb_call() << " calls" << endl;
   cout << "t_global_compute_lambda      = " << t_global_compute_lambda.get_time()      << " s for "<< t_global_compute_lambda.get_nb_call()      << " calls" << endl;
   cout << "t_global_compute_projection  = " << t_global_compute_projection.get_time()  << " s for "<< t_global_compute_projection.get_nb_call()  << " calls" << endl<<endl;
   //cout << "t_global_compute_grad_over_c = " << t_global_compute_grad_over_c.get_time() << " s for "<< t_global_compute_grad_over_c.get_nb_call() << " calls" << endl;
   //cout << "t_global_update_R            = " << t_global_update_R.get_time()            << " s for "<< t_global_update_R.get_nb_call()            << " calls" << endl;
   //cout << "t_global_update_L            = " << t_global_update_L.get_time()            << " s for "<< t_global_update_L.get_nb_call()            << " calls" << endl;
   //cout << "t_check_              = " << t_global_check.get_time()               << " s for "<< t_global_check.get_nb_call()               << " calls" << endl;
   //cout << "t_global_init_fact           = " << t_global_init_fact.get_time()           << " s for "<< t_global_init_fact.get_nb_call()           << " calls" << endl;
   //cout << "t_global_init_fact_from_palm = " << t_global_init_fact_from_palm.get_time() << " s for "<< t_global_init_fact_from_palm.get_nb_call() << " calls" << endl<<endl;
}

void palm4MSA::print_prox_timers() const
{
/*	cout << "prox timers in palm4MSA : " << endl;
   cout << "total t_prox_const  =  " << t_prox_const.get_time()  << " s for "<< nb_call_prox_const  << " calls" << endl;
   cout << "total t_prox_sp  =  " << t_prox_sp.get_time()  << " s for "<< nb_call_prox_sp  << " calls" << endl;
   cout << "total t_prox_spcol  =  " << t_prox_spcol.get_time()  << " s for "<< nb_call_prox_spcol  << " calls" << endl;
   cout << "total t_prox_splin  =  " << t_prox_splin.get_time()  << " s for "<< nb_call_prox_splin  << " calls" << endl;
   cout << "total t_prox_normcol  =  " << t_prox_normcol.get_time()  << " s for "<< nb_call_prox_normcol  << " calls" << endl;	
*/}


void palm4MSA::print_local_timers()const
{
   cout << "timers in palm4MSA : " << endl;
   cout << "t_local_next_step           = " << t_local_next_step.get_time()           << " s for "<< t_local_next_step.get_nb_call()           << " calls" << endl;
   cout << "t grad + updateL + updateR  = " << t_local_update_L.get_time()+t_local_update_R.get_time()+t_local_compute_grad_over_c.get_time()            << " s for "<< t_local_update_L.get_nb_call()            << " calls of grad" << endl;
   cout << "t_local_compute_lambda      = " << t_local_compute_lambda.get_time()      << " s for "<< t_local_compute_lambda.get_nb_call()      << " calls" << endl;
   cout << "t_local_compute_projection  = " << t_local_compute_projection.get_time()  << " s for "<< t_local_compute_projection.get_nb_call()  << " calls" << endl<<endl;
   //cout << "t_local_compute_grad_over_c = " << t_local_compute_grad_over_c.get_time() << " s for "<< t_local_compute_grad_over_c.get_nb_call() << " calls" << endl;
   //cout << "t_local_update_R            = " << t_local_update_R.get_time()            << " s for "<< t_local_update_R.get_nb_call()            << " calls" << endl;
   //cout << "t_local_update_L            = " << t_local_update_L.get_time()            << " s for "<< t_local_update_L.get_nb_call()            << " calls" << endl;
   //cout << "t_check_                    = " << t_local_check.get_time()               << " s for "<< t_local_check.get_nb_call()               << " calls" << endl;
   //cout << "t_local_init_fact           = " << t_local_init_fact.get_time()           << " s for "<< t_local_init_fact.get_nb_call()           << " calls" << endl;
   //cout << "t_local_init_fact_from_palm = " << t_local_init_fact_from_palm.get_time() << " s for "<< t_local_init_fact_from_palm.get_nb_call() << " calls" << endl<<endl;
}







#endif
