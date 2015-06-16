#ifndef __FAUST_PALM4MSA_H__
#define __FAUST_PALM4MSA_H__

#include <iostream>
#include "faust_constant.h"
#include <vector>
#include "faust_mat.h"
#include "stopping_criterion.h"
class faust_constraint_generic;
class faust_params;
class faust_params_palm;

class palm4MSA
{
   public:
      palm4MSA(const faust_params& params_);
      palm4MSA(const faust_params_palm& params_palm_);

      void set_constraint(const std::vector<const faust_constraint_generic*> const_vec_){const_vec=const_vec_;}
      void set_data(const faust_mat& data_){data=data_;}
      void set_nfacts(const int nfact_){nb_fact=nfact_;}
      void set_lambda(faust_real lambda_){lambda = lambda_;}
      void set_lambda(const palm4MSA& palm_){lambda = palm_.lambda;}
      void update_lambda_from_palm(const palm4MSA& palm){lambda *= palm.lambda;}
      faust_real get_lambda()const{return lambda;}
      faust_real get_RMSE()const{return error.norm()/data.getNbRow()/data.getNbCol();}
      const faust_mat& get_res(bool isFactSideLeft_, int ind_)const{return isFactSideLeft_ ? S[0] : S[ind_+1];}
      const faust_mat& get_data()const{return data;}

      void init_fact();      
      void next_step();
      

      ~palm4MSA(){}

   private:
      void check_constraint_validity();
      void init_fact_from_palm(const palm4MSA& palm, bool isFactSideLeft);
      void compute_c();
      void compute_grad_over_c();
      void compute_projection();
      void update_L();
      void update_R();
      void compute_lambda();


   private:
      // RorL_vec matches R if (!isUpdateWayR2L)
      // RorL_vec matches L if (isUpdateWayR2L)
      std::vector<faust_mat> RorL; 

      // LorR_mat matches L if (!isUpdateWayR2L)
      // LorR_mat matches R if (isUpdateWayR2L)
      faust_mat LorR;
      std::vector<faust_mat> S; // contains S_0^i, S_1^i, ...

      
      faust_mat grad_over_c;
      faust_real lipschitz_multiplicator;
      faust_real c; 
      faust_real lambda;
      int ind_fact; //indice de facteur (!= hierarchical_fact::ind_fact : indice de factorisation)
      const bool isUpdateWayR2L;
      const bool verbose;
      faust_mat data;
      faust_mat error; // error = lambda*L*S*R - data
      faust_mat X_hat;
      
      const stopping_criterion stop_crit;

      bool isCComputed;
      bool isGradComputed;
      bool isProjectionComputed;
     
      
      int nb_fact; // number of factors
      std::vector<const faust_constraint_generic*> const_vec; // vector of constraints of size nfact
};

inline void palm4MSA::compute_c()
{ 
   faust_real nL=LorR.spectralNorm();
   faust_real nR=RorL[ind_fact].spectralNorm();
   c=lipschitz_multiplicator*nR*nR*nL*nL;
   std::cout<<"c="<<c<<std::endl;
   isCComputed = true;   
}


#endif
