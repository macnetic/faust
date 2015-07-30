#include "sparse_hierarchical_fact.h"
#include "faust_timer.h"
#include "faust_spmat.h"
#include "faust_core.h"
using namespace std;

//hierarchical_fact::hierarchical_fact(){} // voir avec Luc les parametres par defaut

sparse_hierarchical_fact::sparse_hierarchical_fact(const faust_params& params_):
   ind_fact(0),
   cons(params_.cons),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isFactSideLeft(params_.isFactSideLeft),
   isVerbose(params_.isVerbose),
   nb_fact(params_.nb_fact-1),
   palm_2(sparsePalm4MSA(params_, false)),
   palm_global(sparsePalm4MSA(params_, true)),
   cons_tmp_global(vector<const faust_constraint_generic*>()),
   default_lambda(params_.init_lambda),
   isFactorizationComputed(false),
   errors(std::vector<std::vector<faust_real> >(2,std::vector<faust_real >(params_.nb_fact-1,0.0))){}

void sparse_hierarchical_fact::init()
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
    //cout<<"**************** PALM_GLOBAL ********************"<<endl;
   //palm_global.Display();
   palm_global.init_fact();
    //cout<<"**************** PALM_GLOBAL ********************"<<endl;
   //palm_global.Display();

#ifdef __COMPILE_TIMERS__
t_init.stop();
#endif

}

void sparse_hierarchical_fact::next_step()
{
#ifdef __COMPILE_TIMERS__
t_next_step.start();
#endif

    //cout<<"DEBUT NEXT_STEP"<<endl;  
   if(isFactorizationComputed)
   {
      cerr << "factorization has already been computed" << endl;
      exit(EXIT_FAILURE);
   }

   vector<const faust_constraint_generic*> cons_tmp_2(2);
   cons_tmp_2[0]=cons[0][ind_fact];
   cons_tmp_2[1]=cons[1][ind_fact];
   
   palm_2.set_constraint(cons_tmp_2);
   
   palm_2.init_fact();

   palm_2.set_lambda(default_lambda);
   
   //cout<<"**************** PALM2 ********************"<<endl;
   //palm_2.Display();
   
#ifdef __COMPILE_TIMERS__
palm_2.init_local_timers();
#endif
   //while(palm_2.do_continue())
    //  palm_2.next_step();
 
	palm_2.compute_facts();
//cout<<endl<<endl<<endl; 	
	
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
	//cout<<"**************** PALM_GLOBAL before init_fact ********************"<<endl;
   //palm_global.Display();	

   palm_global.init_fact_from_palm(palm_2, isFactSideLeft);
   //cout<<"**************** PALM_GLOBAL ********************"<<endl;
   //palm_global.Display();	
#ifdef __COMPILE_TIMERS__
palm_global.init_local_timers();
#endif
   //while(palm_global.do_continue())
    //  palm_global.next_step();
	palm_global.compute_facts();
#ifdef __COMPILE_TIMERS__
palm_global.print_local_timers();
#endif
	faust_mat residuum;
	residuum=palm_global.get_residuum(isFactSideLeft); 
    palm_2.set_data(residuum);
   


   //compute_errors();


   ind_fact++;  
    

#ifdef __COMPILE_TIMERS__
t_next_step.stop();
palm_2.print_prox_timers();
#endif
}





void sparse_hierarchical_fact::compute_facts() 
{
   if(isFactorizationComputed)
   {
      cerr << "Error in sparse_hierarchical_fact::compute_facts : factorization has already been computed" << endl;
      exit(EXIT_FAILURE);
   }
  cout<<"compute_fact"<<endl;	
  init();
  for (int i=0 ; i<=nb_fact-1 ; i++)
  {
     cout << "sparse_hierarchical_fact::compute_facts : factorisation "<<i+1<<"/"<<nb_fact <<endl;
     next_step();
  }

  isFactorizationComputed = true;
   
}






#ifdef __COMPILE_TIMERS__
faust_timer sparse_hierarchical_fact::t_init;
faust_timer sparse_hierarchical_fact::t_next_step;

void sparse_hierarchical_fact::print_timers()const
{
   palm_global.print_global_timers();
   cout << "timers in hierarchical_fact :" << endl;
   cout << "t_init      = " << t_init.get_time()      << " s for "<< t_init.get_nb_call()      << " calls" << endl;
   cout << "t_next_step = " << t_next_step.get_time() << " s for "<< t_next_step.get_nb_call() << " calls" << endl<<endl;
}
#endif
