#include "palm4MSA.h"

#include <iostream>
#include "faust_params.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_real.h"
#include "faust_constraint_int.h"
#include "faust_mat.h"
#include "LinAlgebra.h"
#include "prox.h"


palm4MSA::palm4MSA(const faust_params& params_) :
   data(params_.data),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   lambda(params_.init_lambda),
   verbose(params_.isVerbose),
   nb_fact(params_.nb_fact),
   S(params_.init_fact),
   ind_fact(0),
   lipschitz_multiplicator(1.001),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false)
{
   check_constraint_validity();
}

palm4MSA::palm4MSA(const faust_params& params_palm_) :
   data(params_palm_.data),
   isUpdateWayR2L(params_palm_.isUpdateWayR2L),
   lambda(params_palm_.init_lambda),
   verbose(params_palm_.isVerbose),
   stop_crit(params_palm_.stop_crit),
   ind_fact(0),
   lipschitz_multiplicator(1.001),
   isCComputed(false),
   isGradComputed(false),
   isProjectionComputed(false)
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
      faust_mat matrix2project(S[ind_fact]);
      matrix2project -= grad_over_c;
      switch (const_vec[ind_fact]->getConstraintType())
      {
         case CONSTRAINT_NAME_SP:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            prox_sp(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPCOL:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_spcol(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLIN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            prox_splin(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_NORMCOL:
         {
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            prox_normcol(matrix2project, const_real->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLINCOL:
         {
            const faust_constraint_real* const_real = dynamic_cast<const faust_constraint_real*>(const_vec[ind_fact]);
            //prox_splincol(matrix2project, const_real->getParameter());
         }
         break;

         case CONSTRAINT_NAME_L0PEN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l0pen(matrix2project, sqrt(2*const_int->getParameter()/c));
         }
         break;

         case CONSTRAINT_NAME_L1PEN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_l1pen(matrix2project, const_int->getParameter()/c);
         }
         break;

         case CONSTRAINT_NAME_WAV:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_wav(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SP_POS:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp_pos(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_BLKDIAG:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SPLIN_TEST:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_SUPP:
         {
            const faust_constraint_mat* const_mat = dynamic_cast<const faust_constraint_mat*>(const_vec[ind_fact]);
            //prox_sp(matrix2project, const_mat->getParameter());
         }
         break;

         case CONSTRAINT_NAME_NORMLIN:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(matrix2project, const_int->getParameter());
         }
         break;

         case CONSTRAINT_NAME_TOEPLITZ:
         {
            const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
            //prox_sp(matrix2project, const_int->getParameter());
         }
         break;

         default:
            std::cerr << "error in palm4MSA::compute_projection : unknown name of constraint" << std::endl;
            exit(EXIT_FAILURE);
         break;       
      }
   }
   isProjectionComputed = true;
}

void palm4MSA::compute_grad_over_c()
{

   if(!isCComputed) 
   {
      std::cerr << "c must be set before computing grad/c" << std::endl;
      exit(EXIT_FAILURE);
   }

   faust_mat tmp1, tmp2, tmp3, tmp4;
   // tmp1 = L*S
   multiply(L, S[ind_fact], tmp1);

   if (ind_fact == nb_fact-1)
   {
      
      // X_hat = tmp1*R  (= L*S*R )
      multiply(tmp1, R[ind_fact], X_hat);
       
      // error = (X_hat*lambda)-data  (= lambda*L*S*R - data )
      error = X_hat;
      error *= lambda;
      error -= data;
   }
   else
   {
      error = data;
      // error = lambda*tmp1*R - error (= lambda*L*S*R - data )
      gemm(tmp1, R[ind_fact], error, lambda, -1.0, 'N', 'N');
   }
   // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
   gemm(L, error, tmp3, lambda, 0.0, 'T', 'N');
   // grad_over_c = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
   gemm(tmp3, R[ind_fact], grad_over_c,1.0/c, 0.0,'N','T');

   isGradComputed = true;
}

void palm4MSA::compute_lambda()
{
   if (ind_fact != nb_fact-1)
   {
      std::cerr << "error in palm4MSA::compute_lambda : computation of lamda must be done at the end of the iteration through the number of factors" << std::endl;
      exit(EXIT_FAILURE);
   }

   // Xt_Xhat = data'*X_hat
   faust_mat Xt_Xhat;
   gemm(data, X_hat,Xt_Xhat, 1.0, 0.0, 'T','N');
   
   // Xhatt_Xhat = X_hat'*X_hat
   faust_mat Xhatt_Xhat;
   gemm(X_hat, X_hat, Xhatt_Xhat, 1.0, 0.0, 'T','N');

   lambda = Xt_Xhat.trace()/Xhatt_Xhat.trace();
}


void palm4MSA::update_R()
{
   // R[nb_fact-1] est initialise a l'identite lors de la creation de l'objet palm4MSA et n'est pas cense changer
   R[nb_fact-1].setEyes(/*DIMENSION*/);
   for (int i=nb_fact-2 ; i>-1 ; i--)
      //  R[i] = S[i+1] * R[i+1]
      multiply(S[i+1], R[i+1], R[i]);
}


void palm4MSA::check_constraint_validity()
{
   if (nb_fact != S.size())
   {
      std::cerr << "Error in palm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << std::endl;
      exit(EXIT_FAILURE);
   }
}

void palm4MSA::init_fact()
{
   
   if (isUpdateWayR2L)
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

   update_R();
   L.setEyes(/* DIMENSION */);

   for (int j=0 ; j<nb_fact ; j++)
   {
      ind_fact = 0;
      isCComputed = false;
      isGradComputed = false;
      isProjectionComputed = false;
      
      compute_c();
      // X_hat is computed updated by compute_grad_over_c only when j=ind_fact-1
      compute_grad_over_c();
      compute_projection();
      update_L();

      ind_fact++;
   }
   compute_lambda();


}

void palm4MSA::init_fact_from_palm(const palm4MSA& palm2, bool isFactSideLeft)
{
   if (palm2.nb_fact != 2)
   {
      std::cerr << "argument palm2 must contain 2 factors." << std::endl;
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
}

