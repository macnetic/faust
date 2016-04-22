#ifndef __FAUST_PALM4MSA_H__
#define __FAUST_PALM4MSA_H__

#include <iostream>
#include "faust_constant.h"
#include <vector>
#include "stopping_criterion.h"


#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

#include "faust_mat.h"
template<typename T> class faust_mat;
template<typename T> class faust_core;

template<typename T> class faust_constraint_generic;
template<typename T> class faust_params;
template<typename T> class faust_params_palm;
template<typename T> class stopping_criterion;

/*! \class palm4MSA
   * \brief template class implementing palm4MSA (PALM for Multi-layer Sparse Approximation) factorization algorithm
    : <br>
	factorization of a data matrix (dense) into multiple factors (sparse) using PALM
   *
   *\tparam T scalar numeric type, e.g float or double
   */
template<typename T>
class palm4MSA
{

  public:
	/*!
    *  \brief Initialize palm4MSA from faust_params (hierarchical_fact parameter)
	* \tparam params_ : hierarchical_fact parameters
	* \tparam isGlobal_ : if true, the palm4MSA stopping_criterion stop_crit attribute is initialize from params_.stop_crit_global <br> and if false, it is initialize from stop_crit_2facts
    */
    palm4MSA(const faust_params<T>& params_, const bool isGlobal_);
    palm4MSA(const faust_params_palm<T>& params_palm_, const bool isGlobal_=false);

    void set_constraint(const std::vector<const faust_constraint_generic<T>* > const_vec_){const_vec=const_vec_;isConstraintSet=true;}
    void set_data(const faust_mat<T>& data_){data=data_;}
    void set_lambda(const FFPP lambda_){lambda = lambda_;}

	/*! \brief useful in hierarchical_fact, update lambda of palm_global from palm_2 */
    void update_lambda_from_palm(const palm4MSA& palm){lambda *= palm.lambda;}

	/*! \brief Compute the factorization */
    void compute_facts();

	/*! \brief return the multiplicative scalar lambda */
    T get_lambda()const{return lambda;}

    /*! \brief return the RMSE (Root-Mean-Square Error) */
    T get_RMSE()const{return error.norm()/sqrt((double)(data.getNbRow()*data.getNbCol()));}

    const faust_mat<T>& get_res(bool isFactSideLeft_, int ind_)const{return isFactSideLeft_ ? S[0] : S[ind_+1];}
    const faust_mat<T>& get_data()const{return data;}
	void get_facts(faust_core<T> & faust_fact) const;

	/*!
    *  \brief Initialize the factors to the default value, <br>
	* The first factor to be factorized is set to zero matrix whereas all the other are set to identity
    */
    void init_fact(int nb_facts_);
    void next_step();
    /** \brief
     *
     * \param
     * \param
     * \return TRUE or FALSE
     * \warning !!! pre-increment of ind_ite: the value in stop_crit.do_continue is ind_ite+1, not ind_ite
     */
    bool do_continue(){bool cont=stop_crit.do_continue(++ind_ite); if(!cont){ind_ite=-1;isConstraintSet=false;}return cont;}
    //bool do_continue()const{return stop_crit.do_continue(++ind_ite, error);};

	/*!
    *  \brief
	* useful in hierarchical_fact, update the factors of palm_global from palm_2
    */
    void init_fact_from_palm(const palm4MSA& palm, bool isFactSideLeft);
    const std::vector<faust_mat<T> >& get_facts()const {return S;}

    ~palm4MSA(){}

    private:
    void check_constraint_validity();
      void compute_c();
      void compute_grad_over_c();
      void compute_projection();
      void update_L();
      void update_R();
      void compute_lambda();
	  static const char * class_name;
	  static const T lipschitz_multiplicator;

   public:
      stopping_criterion<T> stop_crit;


   private:
      faust_mat<T> data;
      T lambda;
      int nb_fact; // number of factors
      std::vector<faust_mat<T> > S; // contains S_0^i, S_1^i, ...

      // RorL_vec matches R if (!isUpdateWayR2L)
      // RorL_vec matches L if (isUpdateWayR2L)
      std::vector<faust_mat<T> > RorL;
      // LorR_mat matches L if (!isUpdateWayR2L)
      // LorR_mat matches R if (isUpdateWayR2L)
      faust_mat<T> LorR;


      std::vector<const faust_constraint_generic<T>* > const_vec; // vector of constraints of size nfact
      int ind_fact; //indice de facteur (!= hierarchical_fact::ind_fact : indice de factorisation)
      int ind_ite;
      // T lipschitz_multiplicator;
      const bool verbose;
      const bool isUpdateWayR2L;
	  const bool isConstantStepSize;
      bool isCComputed;
      bool isGradComputed;
      bool isProjectionComputed;
      bool isLastFact;
      bool isConstraintSet;
      const bool isGlobal;
      bool isInit; // only used for global factorization (if isGlobal)
      faust_mat<T> grad_over_c;
      T c;
      faust_mat<T> error; // error = lambda*L*S*R - data






#ifdef __COMPILE_TIMERS__
   public:
      faust_timer t_local_compute_projection;
      faust_timer t_local_compute_grad_over_c;
      faust_timer t_local_compute_c;
      faust_timer t_local_compute_lambda;
      faust_timer t_local_update_R;
      faust_timer t_local_update_L;
      faust_timer t_local_check;
      faust_timer t_local_init_fact;
      faust_timer t_local_next_step;
      faust_timer t_local_init_fact_from_palm;


      static faust_timer t_global_compute_projection;
      static faust_timer t_global_compute_grad_over_c;
	  static faust_timer t_global_compute_c;
      static faust_timer t_global_compute_lambda;
      static faust_timer t_global_update_R;
      static faust_timer t_global_update_L;
      static faust_timer t_global_check;
      static faust_timer t_global_init_fact;
      static faust_timer t_global_next_step;
      static faust_timer t_global_init_fact_from_palm;

	  static faust_timer t_prox_const;
	  static faust_timer t_prox_sp;
	  static faust_timer t_prox_spcol;
	  static faust_timer t_prox_splin;
	  static faust_timer t_prox_normcol;

	  static int nb_call_prox_const;
	  static int nb_call_prox_sp;
	  static int nb_call_prox_spcol;
	  static int nb_call_prox_splin;
	  static int nb_call_prox_normcol;







   void init_local_timers();

   void print_global_timers()const;
   void print_local_timers()const;
   void print_prox_timers() const;
#endif

};



#include "palm4MSA.hpp"


#endif
