#ifndef __FAUST_HIERARCHICAL_FACT_H__
#define __FAUST_HIERARCHICAL_FACT_H__

#include "faust_constant.h"
#include <vector>
#include "palm4MSA.h"
#include "faust_params.h"

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

class faust_constraint_generic;


class hierarchical_fact
{
   public:
      //hierarchical_fact(); // voir avec Luc les parametres par defaut
      hierarchical_fact(const faust_params& params_);

      void init();
      void next_step();
      const std::vector<faust_mat>& get_facts()const {return palm_global.get_facts();}
      faust_real get_lambda()const{return palm_global.get_lambda();}


   private:
      const std::vector< std::vector<const faust_constraint_generic*> > cons;
      bool isUpdateWayR2L;
      bool isFactSideLeft; 
      bool isVerbose;
      int ind_fact ; //indice de factorisation (!= palm4MSA::ind_fact : indice de facteur)
      int nb_fact; // nombre de factorisations (!= palm4MSA::nb_fact : nombre de facteurs)
      palm4MSA palm_2;
      palm4MSA palm_global;
      const faust_real default_lambda; // initial value of lambda for factorization into two factors
      std::vector<faust_mat> S;
      std::vector<const faust_constraint_generic*> cons_tmp_global;
     
#ifdef __COMPILE_TIMERS__
   public:
      static faust_timer t_init;
      static faust_timer t_next_step;

    void print_timers()const;
#endif


 
};

#endif
