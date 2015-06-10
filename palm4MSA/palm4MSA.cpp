#include "palm4MSA.h"

#include "faust_params.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_mat.h"
#include "faust_constraint_real.h"
#include "faust_constraint_int.h"
#include "faust_mat.h"

palm4MSA::palm4MSA(const faust_params& params_) :
   data(params_.data),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   lambda(params_.init_lambda),
   verbose(params_.isVerbose),
   ind_fact(0),
   lipschitz_multiplicator(1.001){}


void palm4MSA::compute_projection()
{
   switch (const_vec[ind_fact]->getConstraintType())
   {
      case CONSTRAINT_NAME_SP:
         const faust_constraint_int* const_int = dynamic_cast<const faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_SPCOL:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_spcol(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_SPLIN:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_splin(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_NORMCOL:
         faust_constraint_real* const_real = dynamic_cast<faust_constraint_real*>(const_vec[ind_fact]);
         prox_normcol(S[ind_fact], const_real->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         faust_constraint_real* const_real = dynamic_cast<faust_constraint_real*>(const_vec[ind_fact]);
         prox_splincol(S[ind_fact], const_real->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_L0PEN:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_l0pen(S[ind_fact], sqrt(2*const_int->getParameter()/c), grad/c);
         break;
      case CONSTRAINT_NAME_L1PEN:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_l1pen(S[ind_fact], const_int->getParameter()/c, grad/c);
         break;
      case CONSTRAINT_NAME_CONST:
         faust_constraint_mat* const_mat = dynamic_cast<faust_constraint_mat*>(const_vec[ind_fact]);
         S[ind_fact] = const_mat->getParameter();
         break;
      case CONSTRAINT_NAME_WAV:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_wav(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_SP_POS:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp_pos(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_TYPE_SPLIN_TEST:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_TYPE_SUPP:
         faust_constraint_mat* const_mat = dynamic_cast<faust_constraint_mat*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_mat->getParameter(), grad/c);
         break;
      case CONSTRAINT_TYPE_NORMLIN:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      case CONSTRAINT_TYPE_TOEPLITZ:
         faust_constraint_int* const_int = dynamic_cast<faust_constraint_int*>(const_vec[ind_fact]);
         prox_sp(S[ind_fact], const_int->getParameter(), grad/c);
         break;
      default:
         cerr << "error in palm4MSA::compute_projection : unknown name of constraint" <<endl;
         exit(EXIT_FAILURE);
         break;       
   }
}

void palm4MSA::compute_grad()
{
   faust_mat tmp1, tmp2, tmp3, tmp4;
   // tmp1 = L*S
   faust_multiply(L, S[ind_fact], tmp1);

   error = data;
   // error = lambda*tmp1 - error (= lambda*L*S*R - data )
   faust_multiply('N','N',lambda,-1.0,tmp1, R[ind_fact], error);
   // tmp3 = lambda*L'*error (= lambda*L' * (lambda*L*S*R - data) )
   faust_multiply('T','N',lambda, 0.0,L, error, tmp3);
   // grad = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - data) * R' )
   faust_multiply('N','T',1/c, 0.0,tmp3, R[ind_fact], grad);
}

void palm4MSA::compute_lamda()
{
   if (ind_fact != nb_fact-1)
   {
      std::cerr << "error in palm4MSA::compute_lambda : computation of lamda must be done at the end of the iteration through the number of factors" << endl;
      exit(EXIT_FAILURE);
   }

   faust_mat Xhat;
   faust_multiply(L, S[ind_fact], Xhat);
   set_X_hat(Xhat);

   // Xt_Xhat = data'*Xhat
   faust_mat Xt_Xhat = faust_multiply('T','N',data, Xhat);
   
   // Xhatt_Xhat = Xhat'*Xhat
   Xhatt_Xhat = faust_multiply('T','N',Xhat, Xhat);

   lambda = Xt_Xhat.trace()/Xhatt_Xhat.trace();
}


void palm4MSA::update_R()
{
   // R[nb_facts-1] est initialise a l'identite lors de la creation de l'objet palm4MSA et n'est pas cense changer
   R[nb_facts-1].setIdentity(/*DIMENSION*/);
   for (int i=nb_facts-2 ; i>-1 ; i--)
      //  R[i] = S[i+1] * R[i+1]
      faust_multiply(S[i+1], R[i+1], R[i]);
}


void palm4MSA::check_constraint_validity()
{
   if (nb_fact != S.size())
   {
      std::cerr << "Error in palm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << endl;
      exit(EXIT_FAILURE);
   }
}

void palm4MSA::init_fact()
{
   
   if update_way
      S[0].setZeros(const_vec[0]->getRows(), const_vec[0]->getCols());
      for (int i=1 ; i<nb_fact ; i++)
         S[i].setIdentity(const_vec[i]->getRows(), const_vec[i]->getCols());      
   else
      for (int i=0 ; i<nb_fact-1 ; i++)
         S[i].setIdentity(const_vec[i]->getRows(), const_vec[i]->getCols());      
      S[nb_fact-1].setZeros(const_vec[nb_fact-1]->getRows(), const_vec[nb_fact-1]->getCols());     
}

void palm4MSA::next_step()
{

   update_R();
   L.setIdentity(/* DIMENSION */);

   for (int i=0 ; i<nb_facts ; i++)
   {
      compute_c();
      compute_projection();
      update_L();
   }
   


   ind_fact++;
}

void palm4MSA::init_fact_from_palm(const palm4MSA& palm2, bool isFactSideLeft)
{
   if (palm2.nb_fact != 2)
   {
      std::cerr << "argument palm2 must contain 2 factors." << endl;
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

