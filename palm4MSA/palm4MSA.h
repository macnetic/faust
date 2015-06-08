#ifndef __FAUST_PALM4MSA_H__
#define __FAUST_PALM4MSA_H__

#include "faust_params.h"

class faust_vector;
class faust_matrix;
class faust_spmatrix;
class stopping_criterion;


class palm4MSA
{
   public:
      palm4MSA(const faust_params& params_);

      void compute_grad();
      void compute_lambda();
      inline void compute_c(){ faust_real nL=L.norm(),nR=R[ind_fact].norm();c=lipschitz_multiplicator*nR*nR*nL*nL;}
      void compute_projection();
      void update_R();
      inline void update_L(){L *= S[ind_fact];}
      inline void set_constraint(const vector<faust_constraint>& const_vec_){const_vec=const_vec_;}
      inline void set_data(const faust_mat& data_){data=data_;}
      inline void set_nfacts(int nfact_){nb_fact=nfact_;}
      inline void set_lambda(faust_real lambda_){lambda = lambda_;}
      inline faust_real get_lambda()const{return lambda;}
      inline faust_real get_RMSE()const{return error.norm()/data.getDim1()/data.getDim2();}
      inline const faust_mat& get_res(bool isFactSideLeft_, int ind_)const{return isFactSideLeft_ ? S[0] : S[ind_+1];}
      inline const faust_mat& get_data()const{return data;}
      
      void next_step();
      

      ~palm4MSA(){}

   private:
      void check_constraint_validity();


   private:
      // R : vector containing all posible 
      vecteur<faust_mat> R; 
      faust_mat L;
      vector<faust_mat> S; // contains S_0^i, S_1^i, ...

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
      vector<faust_constraint> const_vec; // vector of constraints of size nfact

}

#endif
