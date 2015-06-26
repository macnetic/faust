#include "palm4MSA_test.h"
#include "faust_mat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_init_from_matio.h"
#include "faust_init_from_matio_mat.h"
#include "palm4MSA.h"
#include "hierarchical_fact.h"
#include "faust_timer.h"
#include <iostream>

using namespace std;


int main()
{
  faust_mat data, init_facts1, init_facts2;

  init_faust_mat_from_matio(data, "config_MEG.mat", "data");
  //init_faust_mat_from_matio(init_facts1, "config_MEG.mat", "init_facts1");
  //init_faust_mat_from_matio(init_facts2, "config_MEG.mat", "init_facts2");

  int cons11_name, cons11_parameter, cons11_row, cons11_col;
  int cons12_name, cons12_parameter, cons12_row, cons12_col;
  int cons13_name, cons13_parameter, cons13_row, cons13_col;
  int cons14_name, cons14_parameter, cons14_row, cons14_col;
  int cons15_name, cons15_parameter, cons15_row, cons15_col;

  int cons21_name, cons21_parameter, cons21_row, cons21_col; 
  int cons22_name, cons22_parameter, cons22_row, cons22_col;
  int cons23_name, cons23_parameter, cons23_row, cons23_col;
  int cons24_name, cons24_parameter, cons24_row, cons24_col;
  int cons25_name, cons25_parameter, cons25_row, cons25_col;


  int nfacts, niter1, niter2;
  bool update_way, verbose, fact_side;
  double init_lambda;

  cons11_name      = init_int_from_matio("config_MEG.mat", "cons11_name");
  cons11_parameter = init_int_from_matio("config_MEG.mat", "cons11_parameter");
  cons11_row       = init_int_from_matio("config_MEG.mat", "cons11_row");
  cons11_col       = init_int_from_matio("config_MEG.mat", "cons11_col");

  cons12_name      = init_int_from_matio("config_MEG.mat", "cons12_name");
  cons12_parameter = init_int_from_matio("config_MEG.mat", "cons12_parameter");
  cons12_row       = init_int_from_matio("config_MEG.mat", "cons12_row");
  cons12_col       = init_int_from_matio("config_MEG.mat", "cons12_col");

  cons13_name      = init_int_from_matio("config_MEG.mat", "cons13_name");
  cons13_parameter = init_int_from_matio("config_MEG.mat", "cons13_parameter");
  cons13_row       = init_int_from_matio("config_MEG.mat", "cons13_row");
  cons13_col       = init_int_from_matio("config_MEG.mat", "cons13_col");

  cons14_name      = init_int_from_matio("config_MEG.mat", "cons14_name");
  cons14_parameter = init_int_from_matio("config_MEG.mat", "cons14_parameter");
  cons14_row       = init_int_from_matio("config_MEG.mat", "cons14_row");
  cons14_col       = init_int_from_matio("config_MEG.mat", "cons14_col");

  cons15_name      = init_int_from_matio("config_MEG.mat", "cons15_name");
  cons15_parameter = init_int_from_matio("config_MEG.mat", "cons15_parameter");
  cons15_row       = init_int_from_matio("config_MEG.mat", "cons15_row");
  cons15_col       = init_int_from_matio("config_MEG.mat", "cons15_col");


  cons21_name      = init_int_from_matio("config_MEG.mat", "cons21_name");
  cons21_parameter = init_double_from_matio("config_MEG.mat", "cons21_parameter");
  cons21_row       = init_int_from_matio("config_MEG.mat", "cons21_row");
  cons21_col       = init_int_from_matio("config_MEG.mat", "cons21_col");

  cons22_name      = init_int_from_matio("config_MEG.mat", "cons22_name");
  cons22_parameter = init_int_from_matio("config_MEG.mat", "cons22_parameter");
  cons22_row       = init_int_from_matio("config_MEG.mat", "cons22_row");
  cons22_col       = init_int_from_matio("config_MEG.mat", "cons22_col");

  cons23_name      = init_int_from_matio("config_MEG.mat", "cons23_name");
  cons23_parameter = init_int_from_matio("config_MEG.mat", "cons23_parameter");
  cons23_row       = init_int_from_matio("config_MEG.mat", "cons23_row");
  cons23_col       = init_int_from_matio("config_MEG.mat", "cons23_col");

  cons24_name      = init_int_from_matio("config_MEG.mat", "cons24_name");
  cons24_parameter = init_int_from_matio("config_MEG.mat", "cons24_parameter");
  cons24_row       = init_int_from_matio("config_MEG.mat", "cons24_row");
  cons24_col       = init_int_from_matio("config_MEG.mat", "cons24_col");

  cons25_name      = init_int_from_matio("config_MEG.mat", "cons25_name");
  cons25_parameter = init_int_from_matio("config_MEG.mat", "cons25_parameter");
  cons25_row       = init_int_from_matio("config_MEG.mat", "cons25_row");
  cons25_col       = init_int_from_matio("config_MEG.mat", "cons25_col");



  nfacts = init_int_from_matio("config_MEG.mat", "nfacts");
  niter1 = init_int_from_matio("config_MEG.mat", "niter1");
  niter2 = init_int_from_matio("config_MEG.mat", "niter2");


  //update_way = init_bool_from_matio("config_MEG.mat", "update_way");
  //verbose = init_bool_from_matio("config_MEG.mat", "verbose");
  //fact_side = init_bool_from_matio("config_MEG.mat", "fact_side");
 
  // Creation du vecteur de contrainte
  const faust_constraint_int  cons11(static_cast<faust_constraint_name>(cons11_name), cons11_parameter, cons11_row, cons11_col);
  const faust_constraint_int  cons12(static_cast<faust_constraint_name>(cons12_name), cons12_parameter, cons12_row, cons12_col);
  const faust_constraint_int  cons13(static_cast<faust_constraint_name>(cons13_name), cons13_parameter, cons13_row, cons13_col);
  const faust_constraint_int  cons14(static_cast<faust_constraint_name>(cons14_name), cons14_parameter, cons14_row, cons14_col);
  const faust_constraint_int  cons15(static_cast<faust_constraint_name>(cons15_name), cons15_parameter, cons15_row, cons15_col);

  const faust_constraint_int  cons21(static_cast<faust_constraint_name>(cons21_name), cons21_parameter, cons21_row, cons21_col);
  const faust_constraint_int  cons22(static_cast<faust_constraint_name>(cons22_name), cons22_parameter, cons22_row, cons22_col);
  const faust_constraint_int  cons23(static_cast<faust_constraint_name>(cons23_name), cons23_parameter, cons23_row, cons23_col);
  const faust_constraint_int  cons24(static_cast<faust_constraint_name>(cons24_name), cons24_parameter, cons24_row, cons24_col);
  const faust_constraint_int  cons25(static_cast<faust_constraint_name>(cons25_name), cons25_parameter, cons25_row, cons25_col);


  vector<const faust_constraint_generic*> cons_tmp1;
  cons_tmp1.push_back(&cons11);
  cons_tmp1.push_back(&cons12);
  cons_tmp1.push_back(&cons13);
  cons_tmp1.push_back(&cons14);
  cons_tmp1.push_back(&cons15);

  vector<const faust_constraint_generic*> cons_tmp2;
  cons_tmp2.push_back(&cons21);
  cons_tmp2.push_back(&cons22);
  cons_tmp2.push_back(&cons23);
  cons_tmp2.push_back(&cons24);
  cons_tmp2.push_back(&cons25);


  vector<vector<const faust_constraint_generic*> >cons;
  cons.push_back(cons_tmp1);
  cons.push_back(cons_tmp2);



  // Creation du critere d'arret
  stopping_criterion crit_2(niter1);
  stopping_criterion crit_global(niter2);

  faust_params params(data, nfacts, cons, vector<faust_mat>(), crit_2, crit_global);

  hierarchical_fact hier_fact(params);

  hier_fact.init();
  //faust_timer t1;
  //t1.start();
  for (int i=0 ; i<=nfacts-2 ; i++)
  {
     //cout<<"i="<<i<<endl;
     hier_fact.next_step();
  }

  //t1.stop();
  //cout <<"total hierarchical fact = "<<t1.get_time()<<endl;
return 0;
}
