#include "hierarchical_fact.h"

hierarchical_fact::hierarchical_fact(){} // voir avec Luc les parametres par defaut

hierarchical_fact::hierarchical_fact(const faust_params& params_):
   ind_fact(0);
   cons(params_.cons),
   isUpdateWayR2L(params_.isUpdateWayR2L),
   isFactSideLeft(params_.isFactSideLeft),
   isVerbose(params_.isVerbose),
   palm_2(palm4MSA(params_)),
   palm_global(palm4MSA(params_)),
   lambda(params_.init_lambda){}

void hierarchical_fact::init()
{
   vector<const faust_constraint_generic*> cons_tmp_2(2,faust_constraint_generic*);
   cons_tmp_2.set_data(cons_tmp_global.get_data());


   vector<const faust_constraint_generic*> cons_tmp_global;
   if(isFactSideLeft)
      cons_tmp_global.push_back(cons[0][ind_fact]);
   else
      cons_tmp_global.push_back(cons[1][ind_fact]);

   


}

void hierarchical_fact::next_step()
{
   cons_tmp_2[0]=cons[0][ind_fact];
   cons_tmp_2[1]=cons[1][ind_fact];
   
   palm_2.set_constraint(cons_tmp_2);
   palm_2.set_data(palm_global.get_res(isFactSideLeft),ind_fact);
   palm_2.init_fact();
   palm_2.set_lambda(1.0);
   

   while(palm_2.stop_crit.do_continue())
      palm_2.next_step();


   
   
   

   palm_global.init_fact_from_palm(palm_2, isFactSideLeft);

   palm_global.set_lambda(palm_2);

   if (isFactSideLeft)
   {
      cons_tmp_global[0]=cons[0][ind_fact];
      vector<const faust_constraint_generic*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+1,cons[1][ind_fact]);
   }
   else
   {
      cons_tmp_global[0]=cons[1][ind_fact];
      vector<const faust_constraint_generic*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+k,cons[0][ind_fact]);      
   }

   palm_global.set_constraint(cons_tmp_global);

   



   while(palm_global.stop_crit.do_continue())
      palm_global.next_step();
   



   ind_fact++;   

}

