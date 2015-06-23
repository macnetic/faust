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

#include <fstream>
#include <iomanip>
#define __SP setprecision(20)<<

using namespace std;

palm4MSA::palm4MSA(const faust_params& params_, const bool isGlobal_) :
   data(params_.data),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   lambda(params_.init_lambda),
   verbose(params_.isVerbose),
   nb_fact(0),
   S(params_.init_fact),
   RorL(vector<faust_mat>(2)),
   ind_fact(0),
   ind_ite(0),
   lipschitz_multiplicator(1.001),
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
   data(params_palm_.data),
   isUpdateWayR2L(params_palm_.isUpdateWayR2L),
   lambda(params_palm_.init_lambda),
   verbose(params_palm_.isVerbose),
   nb_fact(params_palm_.nb_fact),
   S(params_palm_.init_fact),
   RorL(vector<faust_mat>(2)),
   LorR(faust_mat(params_palm_.init_fact[0].getNbRow())),
   stop_crit(params_palm_.stop_crit),
   const_vec(params_palm_.cons),
   ind_fact(0),
   ind_ite(0),
   lipschitz_multiplicator(1.001),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false),
   isLastFact(false),
   isConstraintSet(false),
   isGlobal(isGlobal_)
{
   check_constraint_validity();

}

void palm4MSA::compute_projection()
{
   if (const_vec[ind_fact]->getConstraintType() == CONSTRAINT_NAME_CONST)
   {
         const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
         S[ind_fact] = const_mat->getParameter();
   }
   else
   {
      //faust_mat matrix2project(S[ind_fact]);
      S[ind_fact] -= grad_over_c;
      switch (const_vec[ind_fact]->getConstraintType())
      {
         case CONSTRAINT_NAME_SP:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            prox_sp(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPCOL:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            prox_spcol(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLIN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            prox_splin(S[ind_fact], const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_NORMCOL:
         {
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            prox_normcol(S[ind_fact], const_real->getParameter());
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
            //prox_sp_pos(S[ind_fact], const_int->getParameter());
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
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(S[ind_fact], const_int->getParameter());
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
}

void palm4MSA::compute_grad_over_c()
{

   if(!isCComputed) 
   {
      cerr << "c must be set before computing grad/c" << endl;
      exit(EXIT_FAILURE);
   }

   faust_mat tmp1, tmp2, tmp3, tmp4;
   // tmp1 = L*S
   if (!isUpdateWayR2L)
{
//LorR.print_file("L.dat");
//S[ind_fact].print_file("S.dat");
//tmp1.print_file("tmp1.dat");
      multiply(LorR, S[ind_fact], tmp1);
}
   else
      multiply(RorL[ind_fact], S[ind_fact], tmp1);
   
   error = data;
   // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
   if (!isUpdateWayR2L)
      gemm(tmp1, RorL[ind_fact], error, lambda, -1.0, 'N', 'N');
   else
      gemm(tmp1, LorR, error, lambda, -1.0, 'N', 'N');

   // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
   // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
   if (!isUpdateWayR2L)
   {
      gemm(LorR, error, tmp3, lambda, 0.0, 'T', 'N');
      gemm(tmp3, RorL[ind_fact], grad_over_c,1.0/c, 0.0,'N','T');
   }
   else
   {
      gemm(RorL[ind_fact], error, tmp3, lambda, 0.0, 'T', 'N');
      gemm(tmp3, LorR, grad_over_c,1.0/c, 0.0,'N','T');
   }

   isGradComputed = true;
}

void palm4MSA::compute_lambda()
{
   if (!isLastFact)
   {
      cerr << "error in palm4MSA::compute_lambda : computation of lamda must be done at the end of the iteration through the number of factors" << endl;
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
	cout<<lambda<<endl;
   //cout<<__SP lambda<<endl;
}


void palm4MSA::update_R()
{
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
      LorR.multiplyLeft(S[ind_fact]);
}


void palm4MSA::update_L()
{
   if(!isProjectionComputed){
      cerr << "Projection must be computed before updating L" << endl;
      exit(EXIT_FAILURE);
   }
   if(!isUpdateWayR2L)
      LorR *= S[ind_fact];
   else
   {
      RorL.resize(nb_fact);
      RorL[0].resize(const_vec[0]->getRows());
      RorL[0].setEyes();
      for (int i=0 ; i>nb_fact-2 ; i++)
         //  R[i] = S[i+1] * R[i+1]
         multiply(RorL[i] , S[i], RorL[i+1]);
   }
}

void palm4MSA::check_constraint_validity()
{

   if (nb_fact != S.size())
   {
      cerr << "Error in palm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << endl;
      exit(EXIT_FAILURE);
   }
}

void palm4MSA::init_fact(int nb_facts_)
{
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

    
}

void palm4MSA::next_step()
{
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
      // X_hat is computed by compute_grad_over_c only when j=ind_fact-1
      compute_grad_over_c();
      compute_projection();
   
      if(!isUpdateWayR2L)
         update_L();
      else
         update_R();
    
      if (isLastFact)
         compute_lambda();
   }

   delete[] ind_ptr;
   ind_ptr = NULL;

}

void palm4MSA::init_fact_from_palm(const palm4MSA& palm2, bool isFactSideLeft)
{
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
}

