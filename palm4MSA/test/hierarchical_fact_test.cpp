#include "palm4MSA_test.h"
#include "faust_mat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_init_from_matio.h"
#include "palm4MSA.h"
#include "hierarchical_fact.h"

#include <iostream>

using namespace std;


int main()
{
  faust_mat data, init_facts1, init_facts2;

  init_faust_mat_from_matio_mat(data, "config_compared_palm2.mat", "data");
  init_faust_mat_from_matio_mat(init_facts1, "config_compared_palm2.mat", "init_facts1");
  init_faust_mat_from_matio_mat(init_facts2, "config_compared_palm2.mat", "init_facts2");

  int cons11_name, cons11_parameter, cons11_row, cons11_col;
  int cons12_name, cons12_parameter, cons12_row, cons12_col;
  int cons13_name, cons13_parameter, cons13_row, cons13_col;

faust_real cons21_parameter;
  int cons21_name,                   cons21_row, cons21_col; 
  int cons22_name, cons22_parameter, cons22_row, cons22_col;
  int cons23_name, cons23_parameter, cons23_row, cons23_col;

  faust_real cons2_parameter;
  int nfacts, niter1, niter2;
  bool update_way, verbose, fact_side;
  double init_lambda;

  cons11_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons11_name");
  cons11_parameter = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons11_parameter");
  cons11_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons11_row");
  cons11_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons11_col");

  cons12_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons12_name");
  cons12_parameter = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons12_parameter");
  cons12_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons12_row");
  cons12_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons12_col");

  cons13_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons13_name");
  cons13_parameter = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons13_parameter");
  cons13_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons13_row");
  cons13_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons13_col");

  cons21_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons21_name");
  cons21_parameter = init_faust_mat_from_matio_double("config_compared_hierarchical_fact.mat", "cons21_parameter");
  cons21_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons21_row");
  cons21_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons21_col");

  cons22_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons22_name");
  cons22_parameter = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons22_parameter");
  cons22_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons22_row");
  cons22_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons22_col");

  cons23_name      = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons23_name");
  cons23_parameter = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons23_parameter");
  cons23_row       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons23_row");
  cons23_col       = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "cons23_col");



  nfacts = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "nfacts");
  niter1 = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "niter1");
  niter2 = init_faust_mat_from_matio_int("config_compared_hierarchical_fact.mat", "niter2");


  update_way = init_faust_mat_from_matio_bool("config_compared_hierarchical_fact.mat", "update_way");
  verbose = init_faust_mat_from_matio_bool("config_compared_hierarchical_fact.mat", "verbose");
  fact_side = init_faust_mat_from_matio_bool("config_compared_hierarchical_fact.mat", "fact_side");
 
  // Creation du vecteur de contrainte
  const faust_constraint_int  cons11(static_cast<faust_constraint_name>(cons11_name), cons11_parameter, cons11_row, cons11_col);
  const faust_constraint_int  cons12(static_cast<faust_constraint_name>(cons12_name), cons12_parameter, cons12_row, cons12_col);
  const faust_constraint_int  cons13(static_cast<faust_constraint_name>(cons13_name), cons13_parameter, cons13_row, cons13_col);
  const faust_constraint_real cons21(static_cast<faust_constraint_name>(cons21_name), cons21_parameter, cons21_row, cons21_col);
  const faust_constraint_int  cons22(static_cast<faust_constraint_name>(cons22_name), cons22_parameter, cons22_row, cons22_col);
  const faust_constraint_int  cons23(static_cast<faust_constraint_name>(cons23_name), cons23_parameter, cons23_row, cons23_col);

  vector<const faust_constraint_generic*> cons_tmp1;
  cons_tmp1.push_back(&cons11);
  cons_tmp1.push_back(&cons12);
  cons_tmp1.push_back(&cons13);

  vector<const faust_constraint_generic*> cons_tmp2;
  cons_tmp2.push_back(&cons21);
  cons_tmp2.push_back(&cons22);
  cons_tmp2.push_back(&cons23);

  vector<vector<const faust_constraint_generic*> >cons;
  cons.push_back(cons_tmp1);
  cons.push_back(cons_tmp2);



  // Creation du critere d'arret
  stopping_criterion crit_2(niter1);
  stopping_criterion crit_global(niter2);

  faust_params params(data, nfacts, cons, vector<faust_mat>(), crit_2, crit_global, verbose, update_way, fact_side);

  hierarchical_fact hier_fact(params);

  hier_fact.init();
  for (int i=0 ; i<=nfacts-2 ; i++)
  {
     //cout<<"i="<<i<<endl;
     hier_fact.next_step();
  }


return 0;
}
