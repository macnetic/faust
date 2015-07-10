#include "faust_core.h"
#include "faust_vec.h"
#include "hierarchical_fact.h"
#include "faust_params.h"

using namespace std;



faust_core::faust_core() :
   data(std::vector<faust_spmat>()),
   isDataInit(false),
   isEstimateComputed(false),
   estimate(faust_mat()),
   totalNonZeros(0)
{}

faust_core::faust_core(const std::vector<faust_spmat>& facts, const faust_real lambda_) :
   data(facts),
   lambda(lambda_),
   isDataInit(true),
   estimate(faust_mat()),
   isEstimateComputed(false),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
}

faust_core::faust_core(const faust_params& params) :
   data(std::vector<faust_spmat>()),
   lambda(0.0),
   isDataInit(false),
   estimate(faust_mat()),
   isEstimateComputed(false),
   totalNonZeros(0)
{
   hierarchical_fact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   lambda = hier_fact.get_lambda();

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

const faust_mat& faust_core::get_estimate()
{
   if (!isEstimateComputed)
      compute_estimate();

   return estimate;
}

void faust_core::compute_estimate()
{
   estimate.resize(data[0].getNbRow()); 
   estimate.setEyes();
   for(int i=1 ; i<data.size() ; i++)
      estimate *= data[i];

   estimate *= lambda;
 
   isEstimateComputed = true; 
}

faust_real faust_core::get_lambda()const
{
   if(!isDataInit)
   {
      cerr << "Error in faust_core::get_facts : factors are not available" << endl;
      exit(EXIT_FAILURE);
   }
   return lambda;
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

#if 0
//multiplication de gauche a droite
void faust_core::operator*=(const faust_core&  f)
{

	bool multiplyWayR2L = false;
	
	int nb_factors = data.size() + f.data.size();
	vector<faust_mat_generic*> tmp(nb_factors);


	if (! multiplyWayR2L) //multiplication de gauche a droite
	{
		for (int i=0 ; i<data.size() ; i++)
			tmp[i] = &(data[i]);
		for (int i=0 ; i<f.data.size() ; i++)
			tmp[i+data.size()] = &(f.data[f.data.size()-1-i]);
	}
	else //multiplication de droite a gauche
	{
		for (int i=0 ; i<f.data.size() ; i++)
			tmp[i] = &(f.data[i]);
		for (int i=0 ; i<data.size() ; i++)
			tmp[i+f.data.size()] = &(data[data.size()-1-i]);
	}

	while (tmp.size() > 1)
	{
		for (int j=0 ; j<(tmp.size()-(tmp.size()%2))/2 ; j++)
		{
			if (! multiplyWayR2L) //multiplication de gauche a droite
				(*tmp[j]) = (*tmp[2*j]) * (*tmp[2*j+1]);
			else //multiplication de droite a gauche
				(*tmp[j]) = (*tmp[2*j+1]) * (*tmp[2*j]);
			
			if ((tmp[j]->isSparse) && (!tmp[j]->stillSparse()))
			{
				faust_mat* tmp_mat = new faust_mat(tmp[j]);
				delete tmp[j];
				tmp[j] = tmp_mat;
			}
		}
		tmp.erase(tmp.begin()+(tmp.size()-(tmp.size%2))/2, tmp.begin()+tmp.size()-(tmp.size%2));
	}
	// estimate = *(tmp[0]);
	// isDataInit = true;
}
#endif

