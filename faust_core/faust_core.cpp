#include "faust_core.h"
#include "faust_vec.h"
#include "hierarchical_fact.h"
#include "faust_params.h"

using namespace std;



faust_core::faust_core() :
   data(std::vector<faust_spmat>()),
   isDataInit(false),
   totalNonZeros(0)
{}

faust_core::faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_ /* =1.0 */) :
   data(facts),
   isDataInit(true),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();

   if(lambda_ != 1.0 && data.size()>0)
      (data[0]) *= lambda_;
}

faust_core::faust_core(const faust_params& params) :
   data(std::vector<faust_spmat>()),
   isDataInit(false),
   totalNonZeros(0)
{
   hierarchical_fact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   (data[0]) *= hier_fact.get_lambda();

   isDataInit = true;
}

void faust_core::get_facts(std::vector<faust_spmat>& sparse_facts)const 
{
   if(!isDataInit)
   {
      cerr << "Error in faust_core::get_facts : factors are not available" << endl;
      exit(EXIT_FAILURE);
   }
   sparse_facts = data;
}

faust_mat faust_core::get_product()
{
   if (!isDataInit)
   {
      cerr << "Error in faust_mat faust_core::get_product : data has not been initialized" << endl;
      exit(EXIT_FAILURE);
   }

   faust_mat prod(data[0].getNbRow()); 
   prod.setEyes();
   for(int i=0 ; i<data.size() ; i++)
      prod *= data[i];

   return prod;
}

long long int faust_core::get_total_nnz()const
{
   if(!isDataInit)
   {
      cerr << "Error in faust_core::get_total_nnz : factors are not available" << endl;
      exit(EXIT_FAILURE);
   }

   return totalNonZeros;
}


void faust_core::operator*=(const faust_core&  f)
{
   for (int i=0 ; i<f.size() ; i++)
   {
      data.push_back(f.data[i]);
      totalNonZeros += f.data[i].getNonZeros();
   }
}

void faust_core::operator*=(const faust_spmat& S)
{
   if(S.getNbRow()*S.getNbCol() < 1)
   {
      cerr << "Error in faust_core::operator*= : empty sparse matrix" << endl;
      exit(EXIT_FAILURE);
   }
   data.push_back(S);
   totalNonZeros += S.getNonZeros();
}

