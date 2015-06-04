#include "palm4MSA.h"

palm4MSA::palm4MSA(const faust_params& params_ )
{
}

void palm4MSA::compute_projection()
{
   switch (const_vec[ind_factorization].getConstraintType())
   {
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPCOL:
         prox_spcol(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPLIN:
         prox_splin(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_NORMCOL:
         prox_normcol(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SPLINCOL:
         prox_splincol(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_L0PEN:
         prox_l0pen(S[ind_factorization], sqrt(2*const_vec[ind_factorization].getParameter()/c), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_L1PEN:
         prox_l1pen(S[ind_factorization], const_vec[ind_factorization].getParameter()/c, (1/c)*grad);
         break;
      case CONSTRAINT_NAME_CONST:
         S[ind_factorization] = const_vec[ind_factorization].getParameter();
         break;
      case CONSTRAINT_NAME_WAV:
         prox_wav(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP_POS:
         prox_sp_pos(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_BLKDIAG:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      case CONSTRAINT_NAME_SP:
         prox_sp(S[ind_factorization], const_vec[ind_factorization].getParameter(), (1/c)*grad);
         break;
      default:
         cerr << "error in palm4MSA::compute_projection : unknown name of constraint" <<endl;
         exit(EXIT_FAILURE);
         break;       
   }
}

void palm4MSA::compute_grad(const faust_mat& X)
{
   faust_mat tmp1, tmp2, tmp3, tmp4;
   // tmp1 = L*S
   faust_multiply(L, S[ind_factorization], tmp1);

   tmp2 = X;
   // tmp2 = lambda*tmp1 - tmp2 (= lambda*L*S*R - X )
   faust_multiply('N','N',lambda,-1.0,tmp1, R[ind_factorization], grad_tmp);
   // tmp3 = lambda*L'*tmp2 (= lambda*L' * (lambda*L*S*R - X) )
   faust_multiply('T','N',lambda, 0.0,L, tmp2, tmp3);
   // grad = 1/c*tmp3*R' (= 1/c*lambda*L' * (lambda*L*S*R - X) * R' )
   faust_multiply('N','T',1/c, 0.0,tmp3, R[ind_factorization], grad);
}

void palm4MSA::compute_lamda(const faust_mat& X)
{
   if (ind_factorization != nb_fact-1)
   {
      cerr << "error in palm4MSA::compute_lambda : computation of lamda must be done at the end of the iteration through the number of factors" << endl;
      exit(EXIT_FAILURE);
   }

   faust_mat Xhat;
   faust_multiply(L, S[ind_factorization], Xhat);

   // Xt_X = X
   faust_mat Xt_Xhat = faust_multiply('T','N',X, Xhat);
   
   // Xt_X = X'
   Xhatt_Xhat = faust_multiply('T','N',Xhat, Xhat);

   lamda = Xt_Xhat.trace()/Xhatt_Xhat.trace();
}

void palm4MSA::update_L()
{
   L *= S[ind_factorization];
}

void palm4MSA::update_L()
{
   L *= S[ind_factorization];
}

void palm4MSA::compute_c()
{
   c = lipschitz_multiplicator * norm(R[ind_factorization])*norm(R[ind_factorization]) * norm(L)*norm(L);
}


void palm4MSA::update_R()
{
   // R[nb_facts-1] est initialise a l'identite lors de la creation de l'objet palm4MSA et n'est pas cense changer
   // R[nb_facts-1].setIdentity(/*DIMENSION*/);
   for (int i=nb_facts-2 ; i>-1 ; i--)
      //  R[i] = S[i+1] * R[i+1]
      faust_multiply(S[i+1], R[i+1], R[i]);
}

