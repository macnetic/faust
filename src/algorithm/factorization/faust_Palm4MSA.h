#ifndef __FAUST_PALM4MSA_H__
#define __FAUST_PALM4MSA_H__

#include <iostream>
#include "faust_constant.h"
#include <vector>
#include "faust_StoppingCriterion.h"


#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

#include "faust_MatDense.h"

template<Device DEVICE> class BlasHandle;

namespace Faust
{


    template<typename FPP,Device DEVICE> class MatDense;
    template<typename FPP,Device DEVICE> class Transform;

    template<typename FPP,Device DEVICE> class ConstraintGeneric;
    template<typename FPP,Device DEVICE> class Params;
    template<typename FPP,Device DEVICE> class ParamsPalm;
    template<typename FPP> class StoppingCriterion;


    /*! \class Palm4MSA
       * \brief template class implementing Palm4MSA (PALM for Multi-layer Sparse Approximation) factorization algorithm
        : <br>
        factorization of a data matrix (dense) into multiple factors (sparse) using PALM
       *
       *\tparam FPP scalar numeric type, e.g float or double
       */


    template<typename FPP,Device DEVICE>
    class Palm4MSA
    {

       public:
            /*!
         *  \brief
         * initialize Palm4MSA from Faust::Params (HierarchicalFact parameter)
         *\tparam isGlobal_ : if true, the Palm4MSA StoppingCriterion stop_crit attribute is initialize from params_.stop_crit_global <br> and if false, it is initialize from stop_crit_2facts
         */
          Palm4MSA(const Faust::Params<FPP,DEVICE>& params_, const BlasHandle<DEVICE> blasHandle, const bool isGlobal_);
          Palm4MSA(const Faust::ParamsPalm<FPP,DEVICE>& params_palm_, const BlasHandle<DEVICE> blasHandle, const bool isGlobal_=false);

          void set_constraint(const std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> const_vec_){const_vec=const_vec_;isConstraintSet=true;}
          void set_data(const Faust::MatDense<FPP,DEVICE>& data_){data=data_;}
          void set_lambda(const FPP lambda_){lambda = lambda_;}

          /*!
         *  \brief
         * useful in HierarchicalFact, update lambda of palm_global from palm_2
         */
          void update_lambda_from_palm(const Palm4MSA& palm){lambda *= palm.lambda;}

          /*!
         *  \brief
         * compute the factorisation
         */
          void compute_facts();

          /*!
         *  \brief
         * return the multiplicative scalar lambda
         */
          FPP get_lambda()const{return lambda;}

          FPP get_RMSE()const{return error.norm()/sqrt((double)(data.getNbRow()*data.getNbCol()));}
          const Faust::MatDense<FPP,DEVICE>& get_res(bool isFactSideLeft_, int ind_)const{return isFactSideLeft_ ? S[0] : S[ind_+1];}
          const Faust::MatDense<FPP,DEVICE>& get_data()const{return data;}
          void get_facts(Faust::Transform<FPP,DEVICE> & faust_fact) const;


              /*!
         *  \brief
         * initialize the factors to the default value,
         * the first factor to be factorised is set to zero matrix
         * whereas all the other are set to identity
         */
          void init_fact(int nb_facts_);
          void next_step();
          bool do_continue(){bool cont=stop_crit.do_continue(++ind_ite); if(!cont){ind_ite=-1;isConstraintSet=false;}return cont;} // CAUTION !!! pre-increment of ind_ite: the value in stop_crit.do_continue is ind_ite+1, not ind_ite
          //bool do_continue()const{return stop_crit.do_continue(++ind_ite, error);};

          /*!
         *  \brief
         * useful in HierarchicalFact, update the factors of palm_global from palm_2
         */
          void init_fact_from_palm(const Palm4MSA& palm, bool isFactSideLeft);
          const std::vector<Faust::MatDense<FPP,DEVICE> >& get_facts()const {return S;}

          ~Palm4MSA(){}

       private:
          void check_constraint_validity();
          void compute_c();
          void compute_grad_over_c();
          void compute_projection();
          void update_L();
          void update_R();
          void compute_lambda();
          static const char * class_name;
          static const FPP lipschitz_multiplicator;

       public:
          Faust::StoppingCriterion<FPP> stop_crit;


       private:
          // modif AL AL
          Faust::MatDense<FPP,DEVICE> data;

          FPP lambda;
          int nb_fact; // number of factors
          std::vector<Faust::MatDense<FPP,DEVICE> > S; // contains S_0^i, S_1^i, ...

          // RorL_vec matches R if (!isUpdateWayR2L)
          // RorL_vec matches L if (isUpdateWayR2L)
          std::vector<Faust::MatDense<FPP,DEVICE> > RorL;
          // LorR_mat matches L if (!isUpdateWayR2L)
          // LorR_mat matches R if (isUpdateWayR2L)

          // modif AL AL
          Faust::MatDense<FPP,DEVICE> LorR;
          //Faust::MatDense<FPP,DEVICE>& LorR;



          std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> const_vec; // vector of constraints of size nfact
          int ind_fact; //indice de facteur (!= HierarchicalFact::ind_fact : indice de factorisation)
          int ind_ite;
          // FPP lipschitz_multiplicator;
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
          Faust::MatDense<FPP,DEVICE> grad_over_c;
          FPP c;
          Faust::MatDense<FPP,DEVICE> error; // error = lambda*L*S*R - data
          BlasHandle<DEVICE> blas_handle;






    #ifdef __COMPILE_TIMERS__
       public:
          Faust::Timer_gpu t_local_compute_projection;
          Faust::Timer_gpu t_local_compute_grad_over_c;
          Faust::Timer_gpu t_local_compute_c;
          Faust::Timer_gpu t_local_compute_lambda;
          Faust::Timer_gpu t_local_update_R;
          Faust::Timer_gpu t_local_update_L;
          Faust::Timer_gpu t_local_check;
          Faust::Timer_gpu t_local_init_fact;
          Faust::Timer_gpu t_local_next_step;
          Faust::Timer_gpu t_local_init_fact_from_palm;


          static Faust::Timer_gpu t_global_compute_projection;
          static Faust::Timer_gpu t_global_compute_grad_over_c;
          static Faust::Timer_gpu t_global_compute_c;
          static Faust::Timer_gpu t_global_compute_lambda;
          static Faust::Timer_gpu t_global_update_R;
          static Faust::Timer_gpu t_global_update_L;
          static Faust::Timer_gpu t_global_check;
          static Faust::Timer_gpu t_global_init_fact;
          static Faust::Timer_gpu t_global_next_step;
          static Faust::Timer_gpu t_global_init_fact_from_palm;

          static Faust::Timer_gpu t_prox_const;
          static Faust::Timer_gpu t_prox_sp;
          static Faust::Timer_gpu t_prox_spcol;
          static Faust::Timer_gpu t_prox_splin;
          static Faust::Timer_gpu t_prox_normcol;

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

}

#include "faust_Palm4MSA.hpp"


#endif


