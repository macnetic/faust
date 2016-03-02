#include "faust_cu_mat.h"
#include "faust_cu_spmat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_init_from_matio_params.h"
#include "palm4MSA_cu.h"
#include "hierarchical_fact_cu.h"
#include "faust_cu_timer.h"
#include <iostream>
#include <iomanip>

using namespace std;
typedef double faust_real;//faust floating point precision

int main()
{	
  if (typeid(faust_real) == typeid(double))
  {
	cout<<"floating point precision == double"<<endl;
  }

  if (typeid(faust_real) == typeid(float))
  {
	cout<<"floating point precision == float"<<endl;
  }


   cublasHandle_t cublasHandle;
   cublasStatus_t cublasStat = cublasCreate(&cublasHandle);

   cusparseHandle_t cusparseHandle;
   cusparseStatus_t cusparseStat = cusparseCreate(&cusparseHandle);



  faust_mat<faust_real> data, init_facts1, init_facts2;

  char config_MEG_filename[] = "/home/tgautrai/faust2/test/data/config_MEG.mat";

  
  faust_params<faust_real> params; init_params_from_matiofile(params,config_MEG_filename,"params");
  params.Display();
  cout<<"launch"<<endl;	
  hierarchical_fact_cu<faust_real> hier_fact(params, cublasHandle, cusparseHandle);

  faust_cu_timer t1;
  t1.start();
     
   hier_fact.compute_facts();

  t1.stop();
#ifdef __COMPILE_TIMERS__
  hier_fact.print_timers();
#endif
  cout <<"total hierarchical fact = "<<t1.get_time()<<endl;

  vector<faust_cu_spmat<faust_real> > facts;
  hier_fact.get_facts(facts);
  (facts[0]) *= hier_fact.get_lambda();
  char nomFichier[100];
  for (int i=0 ; i<facts.size() ; i++)
  {
     sprintf(nomFichier, "@FAUST_TESTOUTPUT_BIN_DIR@/facts%d_cpp.dat",i);
     facts[i].print_file(nomFichier);
  }
  cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;

   cusparseDestroy(cusparseHandle);
   cublasDestroy(cublasHandle);

return 0;
}
