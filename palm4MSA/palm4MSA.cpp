#include "palm4MSA.h"

palm4MSA::palm4MSA(const faust_params& params_) :
   data(params_.data),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   lamda(params_.init_lambda),
   verbose(params_.verbose),
   const_vec(const_vec_),
   ind_fact(0),
   lipschitz_multiplicator(1.001){}


void palm4MSA::compute_projection()
{
   switch (const_vec[ind_fact].getConstraintType())
   {
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPCOL:
         prox_spcol(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPLIN:
         prox_splin(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_NORMCOL:
         prox_normcol(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         prox_splincol(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_L0PEN:
         prox_l0pen(S[ind_fact], sqrt(2*const_vec[ind_fact].getParameter()/c), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_L1PEN:
         prox_l1pen(S[ind_fact], const_vec[ind_fact].getParameter()/c, (1/c)*grad);
         break;
      case CONSTRAINT_NAME_CONST:
         S[ind_fact] = const_vec[ind_fact].getParameter();
         break;
      case CONSTRAINT_NAME_WAV:
         prox_wav(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP_POS:
         prox_sp_pos(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_fact], const_vec[ind_fact].getParameter(), (1/c)*grad);
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
      cerr << "error in palm4MSA::compute_lambda : computation of lamda must be done at the end of the iteration through the number of factors" << endl;
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
   // R[nb_facts-1].setIdentity(/*DIMENSION*/);
   for (int i=nb_facts-2 ; i>-1 ; i--)
      //  R[i] = S[i+1] * R[i+1]
      faust_multiply(S[i+1], R[i+1], R[i]);
}


void palm4MSA::check_constraint_validity()
{
   if (nb_fact != S.size())
   {
      cerr << "Error in palm4MSA::check_constraint_validity : Wrong initialization: params.nfacts and params.init_facts are in conflict" << endl;
      exit(EXIT_FAILURE);
   }
}


void palm4MSA::next_step()
{
}
