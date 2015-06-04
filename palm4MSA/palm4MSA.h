#ifndef __palm4MSA_H__
#define __palm4MSA_H__

#include "faust_params.h"

class faust_vector;
class faust_matrix;
class faust_spmatrix;
class stopping_criterion;


class palm4MSA
{
   public:
      palm4MSA(
         const faust_params& params_, 
         vector<faust_constraint> const_vec_,
         ) :
            isUpdateWayR2L(params_.isUpdateWayR2L),
            lamda(params_.init_lambda), 
            const_vec(const_vec_),
            ind_factorization(0),
            lipschitz_multiplicator(1.001){nb_it = ()};

      palm4MSA(const vector<faust_spmat>& S_): palm4MSA(), S(S_);
      palm4MSA(): faust_mat L_, faust_spmat S_;

      void compute_grad(const faust_mat& X);
      void compute_lambda(const faust_mat& X);
      void compute_c(const faust_mat& X );
      void compute_projection();
      void update_R();
      void update_L();
      void compute_X_hat();
      void set_constraint(vector<faust_constraint> const_vec_) {const_vec = const_vec_;}

      ~palm4MSA(){}

   private:
      bool check_constraint_validity();


   private:
      // R : vector containing all posible 
      vecteur<faust_mat> R; 
      faust_mat L;
      vector<faust_spmat> S; // contains S_1^i, S_2^i, ..., S_

      faust_mat grad;
      const faust_real lipschitz_multiplicator;
      faust_real c; 
      faust_real lambda;
      int ind_factorization;
      bool isUpdateWayR2L;
      
      stopping_criterion stop_crit;
      /*const int nb_it;   // number of iterations
      // if isStoppingCriterionError then criterion is error else criterion is number of iteration
      bool  isStoppingCriterionError;
      const faust_real errorThreshold;
      // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than 
      int maxIteration;*/
      
      const int nb_fact; // number of factors
      vector<faust_constraint> const_vec; // vector of constraints of size nfact

}

#endif
