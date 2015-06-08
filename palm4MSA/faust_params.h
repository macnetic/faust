#ifndef __FAUST_PARAMS_H__
#define __FAUST_PARAMS_H__

class faust_mat;
class stopping_criterion;



class faust_params
{
   private:
      // constructeur prive de base. Ajout d'un argument inutile (unused)
      // pour eviter la confusion avec le constructeur qui contient les
      // memes trois premiers parametres plus des parametres par defauts 
      faust_params::faust_params(
         const faust_mat& data_,
         const int nb_fact_,
         const vector<vector<faust_constraint> >& cons_,
         const char unused);

   public:

      faust_params(
         const faust_mat& data_,
         const int nb_fact_,
         const vector<vector<faust_constraint> >& cons_,
         const stopping_criterion& stop_crit_2facts = stopping_criterion() ,
         const stopping_criterion& stop_crit_global = stopping_criterion() ,
         const bool isVerbose_ = false ,
         const bool isUpdateWayR2L_ = false ,
         const bool isFactSideLeft_ = false ,
         const faust_real init_lambda_ = 1.0 ,
         const vector<faust_spmat>& init_fact_,) :
            faust_params(data_, nb_fact_, cons_, '\0'),
            niter1(niter1_),
            niter2(niter2),
            isVerbose(isVerbose_),
            isUpdateWayR2L(isUpdateWayR2L_),
            isFactSideLeft(isFactSideLeft_),
            init_lambda(init_lambda_),
            init_fact(init_fact_) {}

      ~faust_params(){}


   public:
      // Required members
      faust_mat data;
      int nb_fact; // number of factors
      vector<vector<faust_constraint> > cons; // vector of constraints

      // Optional members (set to default values if not defined)
      bool isFactSideLeft;
      bool isVerbose;
      bool isUpdateWayR2L;
      vector<faust_spmat> init_fact;
      faust_real init_lambda;

      const int nb_rows; // number of rows of the first factor
      const int nb_cols; // number of columns of the last factor
     
      stopping_criterion stop_crit_2facts;
      stopping_criterion stop_crit_global;
      /*const int nb_it;   // number of iterations
      // if isStoppingCriterionError then criterion is error else criterion is number of iteration
      bool  isStoppingCriterionError;
      const faust_real errorThreshold;
      // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than 
      int maxIteration;*/
}

#endif
