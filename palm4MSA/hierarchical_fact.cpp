#include "hierarchical_fact.h"
#include "faust_timer.h"
using namespace std;

//hierarchical_fact::hierarchical_fact(){} // voir avec Luc les parametres par defaut

hierarchical_fact::hierarchical_fact(const faust_params& params_):
   ind_fact(0),
   cons(params_.cons),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isFactSideLeft(params_.isFactSideLeft),
   isVerbose(params_.isVerbose),
   nb_fact(params_.nb_fact-1),
   palm_2(palm4MSA(params_, false)),
   palm_global(palm4MSA(params_, true)),
   cons_tmp_global(vector<const faust_constraint_generic*>()),
   default_lambda(params_.init_lambda){}

void hierarchical_fact::init()
{

   cons_tmp_global.clear();
   if(isFactSideLeft)
      cons_tmp_global.push_back(cons[0][ind_fact]);
   else
      cons_tmp_global.push_back(cons[1][ind_fact]);
   

   palm_global.set_constraint(cons_tmp_global);
   palm_global.init_fact(1);
   


}

void hierarchical_fact::next_step()
{
   vector<const faust_constraint_generic*> cons_tmp_2(2);
   cons_tmp_2[0]=cons[0][ind_fact];
   cons_tmp_2[1]=cons[1][ind_fact];
   
   palm_2.set_constraint(cons_tmp_2);
   
   palm_2.init_fact(2);

   palm_2.set_lambda(default_lambda);
   
   //faust_timer t1;
   //t1.start();
   while(palm_2.do_continue())
      palm_2.next_step();
   //t1.stop();
   //cout << "palm2 "<< ind_fact << " = " << t1.get_time()<<endl;

   

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

   //faust_timer t2;
   //t2.start();
   while(palm_global.do_continue())
      palm_global.next_step();
   //t2.stop();
   //cout << "palm_global "<< ind_fact << " = " << t2.get_time()<<endl;

   palm_2.set_data(palm_global.get_res(isFactSideLeft, ind_fact));


   ind_fact++;   

}

