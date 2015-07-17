#include "faust_core.h"
#include "faust_vec.h"
#include "hierarchical_fact.h"
#include "faust_params.h"

using namespace std;



faust_core::faust_core() :
   data(std::vector<faust_spmat>()),
   totalNonZeros(0)
{}

faust_core::faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_ /* =1.0 */) :
   data(facts),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();

   if(lambda_ != 1.0 && data.size()>0)
      (data[0]) *= lambda_;
}

faust_core::faust_core(const faust_params& params) :
   data(std::vector<faust_spmat>()),
   totalNonZeros(0)
{
   hierarchical_fact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   (data[0]) *= hier_fact.get_lambda();

}

faust_mat faust_core::get_product()
{

   faust_mat prod(data[0].getNbRow()); 
   prod.setEyes();
   for(int i=0 ; i<data.size() ; i++)
      prod *= data[i];

   return prod;
}

void faust_core::push_back(const faust_spmat& S)
{
   if (size()>0)
      if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      {
         cerr << "Error in faust_core::push_back : incorrect dimensions" << endl;
         exit(EXIT_FAILURE);
      }
   data.push_back(S);
   totalNonZeros += S.getNonZeros();
	
}

void faust_core::operator*=(const faust_core&  f)
{
   for (int i=0 ; i<f.size() ; i++)
      push_back(f.data[i]);
}


void faust_core::Display()const
{
   

   for (int i=0 ; i<size() ; i++)
   {
      if(data[i].getNbCol()>20)
         cout << "data["<<i<<"] = "<< data[i].getNbRow() <<"x"<< data[i].getNbCol()<<" sparse matrix with "<<data[i].getNonZeros()<<" non-zero values"<<endl;
      else
      {
         cout << "data["<<i<<"] = "<<endl;
         data[i].Display();
         cout<<endl;
      }
   }	
   cout<<endl;
   
}
