#ifndef __FAUST_HIERARCHICAL_FACT_HPP__
#define __FAUST_HIERARCHICAL_FACT_HPP__

//#include "hierarchical_fact.h"
#ifdef __COMPILE_TIMERS__
#include "faust_timer.h"
#endif

#ifdef __COMPILE_GPU__
   #include "faust_cu_spmat.h"
   #include "faust_core_cu.h"
#else
   #include "faust_spmat.h"
   #include "faust_core.h"
#endif

#include "faust_exception.h"
using namespace std;

//hierarchical_fact::hierarchical_fact(){} // voir avec Luc les parametres par defaut

template<typename T>
const char * hierarchical_fact<T>::class_name="hierarchical_fact";

template<typename T>
hierarchical_fact<T>::hierarchical_fact(const faust_params<T>& params_):
   ind_fact(0),
   cons(params_.cons),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isFactSideLeft(params_.isFactSideLeft),
   isVerbose(params_.isVerbose),
   nb_fact(params_.nb_fact-1),
   palm_2(palm4MSA<T>(params_, false)),
   palm_global(palm4MSA<T>(params_, true)),
   cons_tmp_global(vector<const faust_constraint_generic*>()),
   default_lambda(params_.init_lambda),
   isFactorizationComputed(false),
   errors(std::vector<std::vector<T> >(2,std::vector<T >(params_.nb_fact-1,0.0))){}


template<typename T>
void hierarchical_fact<T>::init()
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

template<typename T>
void hierarchical_fact<T>::next_step()
{
#ifdef __COMPILE_TIMERS__
t_next_step.start();
#endif


    if(isFactorizationComputed)
    {
        handleError(class_name,"next_step : factorization has already been computed");
    }

    vector<const faust_constraint_generic*> cons_tmp_2(2);
    cons_tmp_2[0]=cons[0][ind_fact];
    cons_tmp_2[1]=cons[1][ind_fact];


    palm_2.set_constraint(cons_tmp_2); /**< setting the constraint parameters*/

    palm_2.init_fact(2);
    palm_2.set_lambda(default_lambda);

#ifdef __COMPILE_TIMERS__
palm_2.init_local_timers();
#endif
   //while(palm_2.do_continue())
    //  palm_2.next_step();
	palm_2.compute_facts(); /**< Compute the hierarchical factorization of dense matrix. */


#ifdef __COMPILE_TIMERS__
palm_2.print_local_timers();
#endif
    palm_global.update_lambda_from_palm(palm_2);


    if (isFactSideLeft)
    {
        cons_tmp_global[0]=cons[0][ind_fact];
        vector<const faust_constraint_generic*>::iterator it;
        it = cons_tmp_global.begin();
        cons_tmp_global.insert(it+1,cons[1][ind_fact]);
    }
    else
    {
        vector<const faust_constraint_generic*>::iterator it;
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

template<typename T>
void hierarchical_fact<T>::get_facts(faust_core<T> & fact)const
{
	std::vector<faust_spmat<T> > spfacts;
	get_facts(spfacts);
	faust_core<T> res(spfacts);
	fact = res;
}

template<typename T>
void hierarchical_fact<T>::get_facts(std::vector<faust_spmat<T> >& sparse_facts)const
{
   /*if(!isFactorizationComputed)
   {
      cerr << "Error in hierarchical_fact<T>::get_facts : factorization has not been computed" << endl;
      exit(EXIT_FAILURE);
   }*/

   const std::vector<faust_mat<T> >& full_facts = palm_global.get_facts();
   sparse_facts.resize(full_facts.size());
   for (int i=0 ; i<sparse_facts.size() ; i++)
      sparse_facts[i] = full_facts[i];
}



template<typename T>
void hierarchical_fact<T>::compute_facts()
{
   if(isFactorizationComputed)
   {
      handleError(class_name,"compute_facts : factorization has already been computed");
   }

  init();
  for (int i=0 ; i<=nb_fact-1 ; i++)
  {
     cout << "hierarchical_fact<T>::compute_facts : factorisation "<<i+1<<"/"<<nb_fact <<endl;
     next_step();
  }

  isFactorizationComputed = true;

}


template<typename T>
const std::vector<std::vector< T> >& hierarchical_fact<T>::get_errors()const
{
    if(!isFactorizationComputed)
    {
        handleError(class_name,"get_errors() : Factorization has not been computed");
    }
    return errors;
}

template<typename T>
void hierarchical_fact<T>::compute_errors()
{
    vector<faust_spmat<T> > sp_facts;
    get_facts(sp_facts);


    faust_core<T> faust_core_tmp(sp_facts, get_lambda());
    const faust_mat<T> estimate_mat = faust_core_tmp.get_product();

    faust_mat<T> data(palm_global.get_data());

    T data_norm = data.norm();
    data -= estimate_mat;
    errors[0][ind_fact] =  estimate_mat.norm()/data_norm;
    errors[1][ind_fact] =  faust_core_tmp.get_total_nnz()/data.getNbRow()/data.getNbCol();

}


#ifdef __COMPILE_TIMERS__
template<typename T> faust_timer hierarchical_fact<T>::t_init;
template<typename T> faust_timer hierarchical_fact<T>::t_next_step;

template<typename T>
void hierarchical_fact<T>::print_timers()const
{
   palm_global.print_global_timers();
   cout << "timers in hierarchical_fact :" << endl;
   cout << "t_init      = " << t_init.get_time()      << " s for "<< t_init.get_nb_call()      << " calls" << endl;
   cout << "t_next_step = " << t_next_step.get_time() << " s for "<< t_next_step.get_nb_call() << " calls" << endl<<endl;
}
#endif

#endif
