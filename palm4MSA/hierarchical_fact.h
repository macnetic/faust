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
template<typename T> class palm4MSA;
template<typename T> class faust_mat;
template<typename T> class faust_spmat;
template<typename T> class faust_core;
template<typename T> class stopping_criterion;


template<typename T>
class hierarchical_fact
{
   public:
      //hierarchical_fact(); // voir avec Luc les parametres par defaut
      hierarchical_fact(const faust_params<T>& params_);

      void get_facts(std::vector<faust_spmat<T> >&)const;
	  void get_facts(std::vector<faust_mat<T> >& fact)const{fact = palm_global.get_facts();}
      void compute_facts();
      T get_lambda()const{return palm_global.get_lambda();}
      const std::vector<std::vector< T> >& get_errors()const;


private:
      void init();
      void next_step();
      void compute_errors();


   private:
      const std::vector< std::vector<const faust_constraint_generic*> > cons;
      bool isUpdateWayR2L;
      bool isFactSideLeft; 
      bool isVerbose;
      int ind_fact ; //indice de factorisation (!= palm4MSA::ind_fact : indice de facteur)
      int nb_fact; // nombre de factorisations (!= palm4MSA::nb_fact : nombre de facteurs)
      palm4MSA<T> palm_2;
      palm4MSA<T> palm_global;
      const T default_lambda; // initial value of lambda for factorization into two factors
      //std::vector<faust_mat<T> > S;
      std::vector<const faust_constraint_generic*> cons_tmp_global;
      bool isFactorizationComputed;
      std::vector<std::vector<T> > errors;
	  static const char * class_name;
      
     
#ifdef __COMPILE_TIMERS__
   public:
      static faust_timer t_init;
      static faust_timer t_next_step;

    void print_timers()const;
	//void print_prox_timers()const;
#endif


 
};


#include "hierarchical_fact.hpp"

#endif
