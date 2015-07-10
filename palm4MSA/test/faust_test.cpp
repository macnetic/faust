#include "faust_core.h"
#include "faust_params.h"
#include "faust_init_from_matio_params.h"
#include "hierarchical_fact.h"
#include "faust_vec.h"
#include "faust_timer.h"

#include <iostream>

using namespace std;


int main()
{

bool testMEG  = false;


  // init from params

  faust_params params;
  if (testMEG)
     init_params_from_matiofile(params, "config_MEG.mat", "params");
  else
     init_params_from_matiofile(params, "config_compared_hierarchical_fact.mat", "params");
    
//cout << "update_way=" << params.isUpdateWayR2L << endl;
 
  vector<faust_spmat> facts;
  faust_core faust(params);
  faust.get_facts(facts);
  faust_real lambda = faust.get_lambda(); 

cout << "lambda=" << lambda << endl;

  char filename[100];
  for (int i=0 ; i<facts.size() ; i++)
  {
     if (testMEG)
        sprintf(filename, "facts_MEG%d.dat", i+1);
     else
        sprintf(filename, "facts_hier_fact%d.dat", i+1);
     facts[i].print_file(filename);
  }

  // init from vector of faust_spmat

  /*vector<faust_spmat> facts(6);
  char filename[100];
  for (int i=0 ; i<facts.size() ; i++)
  {
     if (testMEG)
        sprintf(filename, "facts_MEG%d.txt", i);
     else
        sprintf(filename, "facts_hier_fact%d.txt", i);
     facts[i].init_from_file(filename);
  }
  faust_core faust(facts, lambda);*/




  faust_vec vec_in(facts[facts.size()-1].getNbCol());
  for (int i=0 ; i<vec_in.size() ; i++ )
     vec_in[i] = i*0.015 ;
     
  faust_vec vec_out_faust, vec_out_mat;
  faust_mat full(params.data);

  faust_timer t1;
  t1.start();
     vec_out_faust = faust * vec_in;
  t1.stop();
  cout << "temps multiplication faust-vector = " << t1.get_time() << "s." << endl;

  faust_timer t2;

  t2.start();
     vec_in.multiplyLeft(full);
  t2.stop();
  cout << "temps multiplication matrix-vector = " << t2.get_time() << "s." << endl;
  vec_out_mat = vec_in;

  if (testMEG)
  {
     vec_out_faust.print_file("vec_out_MEG_faust.dat");
     vec_out_mat.print_file("vec_out_MEG_mat.dat");
  }
  else
  {
     vec_out_faust.print_file("vec_out_hier_fact_faust.dat");
     vec_out_mat.print_file("vec_out_hier_fact_mat.dat");
  }

   
 
}


