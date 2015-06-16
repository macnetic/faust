#include "palm4MSA_test.h"
#include "faust_mat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_init_from_matio.h"
#include "palm4MSA.h"

#include <iostream>

using namespace std;


int main()
{
  faust_mat data, init_facts1, init_facts2;

  init_faust_mat_from_matio_mat(data, "config_compared_palm2.mat", "data");
  init_faust_mat_from_matio_mat(init_facts1, "config_compared_palm2.mat", "init_facts1");
  init_faust_mat_from_matio_mat(init_facts2, "config_compared_palm2.mat", "init_facts2");

  int cons1_name, cons1_parameter, cons1_row, cons1_col;
  int cons2_name, cons2_row, cons2_col;
  faust_real cons2_parameter;
  int nfacts, niter;
  bool update_way, verbose;
  double init_lambda;

  cons1_name = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons1_name");
  cons1_parameter = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons1_parameter");
  cons1_row = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons1_row");
  cons1_col = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons1_col");

  cons2_name = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons2_name");
  cons2_parameter = (faust_real) init_faust_mat_from_matio_double("config_compared_palm2.mat", "cons2_parameter");
  cons2_row = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons2_row");
  cons2_col = init_faust_mat_from_matio_int("config_compared_palm2.mat", "cons2_col");


  init_lambda = init_faust_mat_from_matio_double("config_compared_palm2.mat", "init_lambda");
  nfacts = init_faust_mat_from_matio_int("config_compared_palm2.mat", "nfacts");
  niter = init_faust_mat_from_matio_int("config_compared_palm2.mat", "niter");

  update_way = init_faust_mat_from_matio_bool("config_compared_palm2.mat", "update_way");
  verbose = init_faust_mat_from_matio_bool("config_compared_palm2.mat", "verbose");
 
  // Creation du vecteur de contrainte
  const faust_constraint_int cons1(static_cast<faust_constraint_name>(cons1_name), cons1_parameter, cons1_row, cons1_col);
  const faust_constraint_real cons2(static_cast<faust_constraint_name>(cons2_name), cons2_parameter, cons2_row, cons2_col);

  vector<const faust_constraint_generic*> cons;
  cons.push_back(&cons1);
  cons.push_back(&cons2);



  // Creation du vecteur de matrice init_fact;
  vector<faust_mat> init_fact;
  init_fact.push_back(init_facts1);
  init_fact.push_back(init_facts2);

  // Creation du critere d'arret
  stopping_criterion crit(niter);

  faust_params_palm params(data, nfacts, cons, init_fact, crit, verbose, update_way, init_lambda);

  palm4MSA palm2(params);

  palm2.next_step();






return 0;
}
