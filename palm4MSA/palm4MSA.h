#ifndef __FAUST_PALM4MSA_H__
#define __FAUST_PALM4MSA_H__

#include "faust_mat.h"
#include "faust_constraint.h"

class faust_params;
class faust_real;
class stopping_criterion;

class palm4MSA
{
   public:
      palm4MSA(const faust_params& params_);

      void compute_grad();
      void compute_lambda();
      void compute_c(){ faust_real nL=L.norm(),nR=R[ind_fact].norm();c=lipschitz_multiplicator*nR*nR*nL*nL;}
      void compute_projection();
      void update_R();
      void update_L(){L *= S[ind_fact];}
      void set_constraint(const std::vector<const faust_constraint_generic*> const_vec_){const_vec=const_vec_;}
      void set_data(const faust_mat& data_){data=data_;}
      void set_nfacts(const int nfact_){nb_fact=nfact_;}
      void set_lambda(faust_real lambda_){lambda = lambda_;}
      void update_lambda_from_palm(const palm4MSA& palm){lambda *= palm.lambda;}
      faust_real get_lambda()const{return lambda;}
      faust_real get_RMSE()const{return error.norm()/data.getDim1()/data.getDim2();}
      const faust_mat& get_res(bool isFactSideLeft_, int ind_)const{return isFactSideLeft_ ? S[0] : S[ind_+1];}
      const faust_mat& get_data()const{return data;}

      void init_fact();      
      void next_step();
      

      ~palm4MSA(){}

   private:
      void check_constraint_validity();
      void init_fact_from_palm(const palm4MSA& palm, bool isFactSideLeft);


   private:
      // R : vector containing all posible 
      std::vector<faust_mat> R; 
      faust_mat L;
      std::vector<faust_mat> S; // contains S_0^i, S_1^i, ...

      faust_mat grad;
      faust_real lipschitz_multiplicator;
      faust_real c; 
      faust_real lambda;
      int ind_fact; //indice de facteur (!= hierarchical_fact::ind_fact : indice de factorisation)
      const bool isUpdateWayR2L;
      const bool verbose;
      faust_mat data;
      faust_mat error; // error = lambda*L*S*R - data
      
      const stopping_criterion stop_crit;
     
      
      int nb_fact; // number of factors
      std::vector<const faust_constraint_generic*> const_vec; // vector of constraints of size nfact

}

#endif
