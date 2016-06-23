#ifndef __FAUST_HIERARCHICAL_FACT_H__
#define __FAUST_HIERARCHICAL_FACT_H__

#include "faust_constant.h"
#include <vector>
#include "faust_Palm4MSA.h"
#include "faust_Params.h"

#ifdef __COMPILE_FPPIMERS__
  #include "faust_Timer.h"
#endif

/*! \class HierarchicalFact
   * \brief template class implementing hierarchical factorisation algorithm
    : <br>  factorization of a data matrix into multiple factors using Faust::Palm4MSA algorithm mixed with a hierarchical approach.
   *\tparam FPP scalar numeric type, e.g float or double
   */


namespace Faust
{

    template<typename FPP,Device DEVICE> class ConstraintGeneric;
    template<typename FPP,Device DEVICE> class Palm4MSA;
    template<typename FPP,Device DEVICE> class MatDense;
    template<typename FPP,Device DEVICE> class MatSparse;
    template<typename FPP,Device DEVICE> class Transform;

    template<typename FPP> class StoppingCriterion;
//    template<Device DEVICE> class BlasHandle;
//    template<Device DEVICE> class SpBlasHandle;

    template<typename FPP,Device DEVICE>
    class HierarchicalFact
    {

       public:

          HierarchicalFact(const Faust::Params<FPP,DEVICE>& params_, Faust::BlasHandle<DEVICE> cublasHandle, SpBlasHandle<DEVICE> cusparseHandle);
          void get_facts(Faust::Transform<FPP,DEVICE> &)const;
          void get_facts(std::vector<Faust::MatSparse<FPP,DEVICE> >&)const;
          void get_facts(std::vector<Faust::MatDense<FPP,DEVICE> >& fact)const{fact = palm_global.get_facts();}
          void compute_facts();
          FPP get_lambda()const{return palm_global.get_lambda();}
          const std::vector<std::vector< FPP> >& get_errors()const;


        private:
          void init();
          void next_step();
          void compute_errors();


        private:
          const std::vector< std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> > cons;
          bool m_isUpdateWayR2L;
          bool m_isFactSideLeft;
          bool m_isVerbose;
          int m_indFact ; //indice de factorisation (!= Faust::Palm4MSA::m_indFact : indice de facteur)
          int nbFact; // nombre de factorisations (!= Faust::Palm4MSA::nbFact : nombre de facteurs)
          Faust::Palm4MSA<FPP,DEVICE> palm_2;
          Faust::Palm4MSA<FPP,DEVICE> palm_global;
          const FPP default_lambda; // initial value of lambda for factorization into two factors
          //std::vector<Faust::MatDense<FPP,DEVICE> > S;
          std::vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> cons_tmp_global;
          bool isFactorizationComputed;
          std::vector<std::vector<FPP> > errors;
          static const char * m_className;
          Faust::BlasHandle<DEVICE> cublas_handle;
          Faust::SpBlasHandle<DEVICE> cusparse_handle;


    #ifdef __COMPILE_TIMERS__
       public:
          static Faust::Timer t_init;
          static Faust::Timer t_next_step;

        void print_timers()const;
        //void print_prox_timers()const;
    #endif



    };

}

#include "faust_HierarchicalFact.hpp"

#endif
