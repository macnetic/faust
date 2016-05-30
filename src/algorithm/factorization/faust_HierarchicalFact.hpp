#ifndef __HIERARCHICAL_FACT_CU_HPP__
#define __HIERARCHICAL_FACT_CU_HPP__

//#include "faust_HierarchicalFact.h"
#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

#ifdef __COMPILE_GPU__
	#include "faust_MatSparse_gpu.h"
	#include "faust_Transform_gpu.h"
#else
	#include "faust_MatSparse.h"
	#include "faust_Transform.h"
#endif

#include "faust_exception.h"
using namespace std;

//Faust::HierarchicalFact::Faust::HierarchicalFact(){} // voir avec Luc les parametres par defaut

template<typename FPP,Device DEVICE>
const char * Faust::HierarchicalFact<FPP,DEVICE>::class_name="Faust::HierarchicalFact";

template<typename FPP,Device DEVICE>
Faust::HierarchicalFact<FPP,DEVICE>::HierarchicalFact(const Faust::Params<FPP,DEVICE>& params_, BlasHandle<DEVICE> cublasHandle, SpBlasHandle<DEVICE> cusparseHandle):
   ind_fact(0),
   cons(params_.cons),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isFactSideLeft(params_.isFactSideLeft),
   isVerbose(params_.isVerbose),
   nb_fact(params_.nb_fact-1),
   palm_2(Palm4MSA<FPP,DEVICE>(params_, cublasHandle, false)),
   palm_global(Palm4MSA<FPP,DEVICE>(params_, cublasHandle, true)),
   cons_tmp_global(vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>()),
   default_lambda(params_.init_lambda),
   isFactorizationComputed(false),
   errors(std::vector<std::vector<FPP> >(2,std::vector<FPP >(params_.nb_fact-1,0.0))),
   cublas_handle(cublasHandle),
   cusparse_handle(cusparseHandle){}


template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::init()
{
#ifdef __COMPILE_TIMERS__
t_init.start();
#endif

   cons_tmp_global.clear();
   if(isFactSideLeft)
      cons_tmp_global.push_back(cons[0][ind_fact]);
   else
      cons_tmp_global.push_back(cons[1][ind_fact]);


   palm_global.set_constraint(cons_tmp_global);
   palm_global.init_fact(1);


#ifdef __COMPILE_TIMERS__
t_init.stop();
#endif

}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::next_step()
{
#ifdef __COMPILE_TIMERS__
t_next_step.start();
#endif

   if(isFactorizationComputed)
   {
      handleError(class_name,"next_step : factorization has already been computed");
   }

   vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> cons_tmp_2(2);
   cons_tmp_2[0]=cons[0][ind_fact];
   cons_tmp_2[1]=cons[1][ind_fact];


   palm_2.set_constraint(cons_tmp_2);

   palm_2.init_fact(2);

   palm_2.set_lambda(default_lambda);

#ifdef __COMPILE_TIMERS__
palm_2.init_local_timers();
#endif
   //while(palm_2.do_continue())
    //  palm_2.next_step();
	palm_2.compute_facts();


#ifdef __COMPILE_TIMERS__
palm_2.print_local_timers();
#endif
   palm_global.update_lambda_from_palm(palm_2);


   if (isFactSideLeft)
   {
      cons_tmp_global[0]=cons[0][ind_fact];
      typename vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+1,cons[1][ind_fact]);
   }
   else
   {
      typename vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+ind_fact,cons[0][ind_fact]);
      cons_tmp_global[ind_fact+1]=cons[1][ind_fact];
   }

   palm_global.set_constraint(cons_tmp_global);


   palm_global.init_fact_from_palm(palm_2, isFactSideLeft);

#ifdef __COMPILE_TIMERS__
palm_global.init_local_timers();
#endif
   //while(palm_global.do_continue())
    //  palm_global.next_step();
	palm_global.compute_facts();
#ifdef __COMPILE_TIMERS__
palm_global.print_local_timers();
#endif

   palm_2.set_data(palm_global.get_res(isFactSideLeft, ind_fact));


   compute_errors();


   ind_fact++;


#ifdef __COMPILE_TIMERS__
t_next_step.stop();
palm_2.print_prox_timers();
#endif
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::get_facts(Faust::Transform<FPP,DEVICE> & fact)const
{
	std::vector<Faust::MatSparse<FPP,DEVICE> > spfacts;
	get_facts(spfacts);
	Faust::Transform<FPP,DEVICE> res(spfacts);
	fact = res;
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::get_facts(std::vector<Faust::MatSparse<FPP,DEVICE> >& sparse_facts)const
{
   /*if(!isFactorizationComputed)
   {
      cerr << "Error in Faust::HierarchicalFact<FPP,DEVICE>::get_facts : factorization has not been computed" << endl;
      exit(EXIT_FAILURE);
   }*/

   const std::vector<Faust::MatDense<FPP,DEVICE> >& full_facts = palm_global.get_facts();
   sparse_facts.resize(full_facts.size());
   for (int i=0 ; i<sparse_facts.size() ; i++)
{
      sparse_facts[i].init(full_facts[i],cusparse_handle);
}
}



template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::compute_facts()
{
   if(isFactorizationComputed)
   {
      handleError(class_name,"compute_facts : factorization has already been computed");
   }

  init();
  for (int i=0 ; i<=nb_fact-1 ; i++)
  {
     cout << "Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation "<<i+1<<"/"<<nb_fact <<endl;
     next_step();
  }

  isFactorizationComputed = true;

}


template<typename FPP,Device DEVICE>
const std::vector<std::vector< FPP> >& Faust::HierarchicalFact<FPP,DEVICE>::get_errors()const
{
    if(!isFactorizationComputed)
    {
        handleError(class_name,"get_errors() : Factorization has not been computed");
    }
    return errors;
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::compute_errors()
{
   vector<Faust::MatSparse<FPP,DEVICE> > sp_facts;
   get_facts(sp_facts);



   Faust::Transform<FPP,DEVICE> faust_Transform_tmp(sp_facts, get_lambda());
   const Faust::MatDense<FPP,DEVICE> estimate_mat = faust_Transform_tmp.get_product(cublas_handle, cusparse_handle);

   Faust::MatDense<FPP,DEVICE> data(palm_global.get_data());

   FPP data_norm = data.norm();


   data -= estimate_mat;

   errors[0][ind_fact] =  estimate_mat.norm()/data_norm;
   errors[1][ind_fact] =  faust_Transform_tmp.get_total_nnz()/data.getNbRow()/data.getNbCol();

}


#ifdef __COMPILE_TIMERS__
template<typename FPP,Device DEVICE> Faust::Timer Faust::HierarchicalFact<FPP,DEVICE>::t_init;
template<typename FPP,Device DEVICE> Faust::Timer Faust::HierarchicalFact<FPP,DEVICE>::t_next_step;

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::print_timers()const
{
   palm_global.print_global_timers();
   cout << "timers in Faust::HierarchicalFact :" << endl;
   cout << "t_init      = " << t_init.get_time()      << " s for "<< t_init.get_nb_call()      << " calls" << endl;
   cout << "t_next_step = " << t_next_step.get_time() << " s for "<< t_next_step.get_nb_call() << " calls" << endl<<endl;

Faust::MatDense<FPP,DEVICE> M_tmp;
Faust::MatSparse<FPP,DEVICE> S_tmp;
M_tmp.print_timers();
S_tmp.print_timers();

}
#endif

#endif
