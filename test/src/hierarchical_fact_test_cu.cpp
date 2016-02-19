#include "faust_cu_mat.h"
#include "faust_cu_spmat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_real.h"
#include "faust_params.h"
#include "stopping_criterion.h"
#include "hierarchical_fact_cu.h"
#include "faust_timer.h"
#include "faust_init_from_matio_params.h"
#include <iostream>
#include "faust_constant.h"
#include <iomanip>
#include "cublas_v2.h"
#include "cusparse.h"


using namespace std;
typedef float FPP;// floating point precision

int main()
{

  if (typeid(FPP) == typeid(double))
  {
	cout<<"floating point precision == double"<<endl;
  }
  
  if (typeid(FPP) == typeid(float))
  {
	cout<<"floating point precision == float"<<endl;
  }
		
  
   cublasHandle_t cublasHandle;
   cublasStatus_t cublasStat = cublasCreate(&cublasHandle);

   cusparseHandle_t cusparseHandle;
   cusparseStatus_t cusparseStat = cusparseCreate(&cusparseHandle);

                                 
  
	
    // char config_filename[] = "@FAUST_TESTDATA_SRC_DIR@/config_compared_hierarchical_fact.mat";
	 char config_filename[] = "../data/config_compared_hierarchical_fact.mat";

  faust_params<FPP> params;
  init_params_from_matiofile(params,config_filename,"params");
  params.Display();	
  hierarchical_fact_cu<FPP> hier_fact(params, cublasHandle, cusparseHandle);

  faust_timer t1;
  t1.start();

  hier_fact.compute_facts();

  t1.stop();
#ifdef __COMPILE_TIMERS__
  hier_fact.print_timers();
  //hier_fact.print_prox_timers();
#endif
  cout <<"total hierarchical fact = "<<t1.get_time()<<endl;

  vector<faust_cu_spmat<FPP> > facts;
  hier_fact.get_facts(facts);
  (facts[0]) *= hier_fact.get_lambda();
  char nomFichier[100];
  for (unsigned int i=0 ; i<facts.size() ; i++)
  {
     //sprintf(nomFichier, "@FAUST_TESTOUTPUT_BIN_DIR@/facts%d_cpp.dat",i);
     sprintf(nomFichier, "hier_fact_cu_facts%d_cpp.dat",i);
     facts[i].print_file(nomFichier);
  }
  cout<<"lambda="<<std::setprecision(20)<<hier_fact.get_lambda()<<endl;




   cusparseDestroy(cusparseHandle);
   cublasDestroy(cublasHandle);

return 0;
}
