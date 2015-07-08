#include "faust_core.h"
#include "faust_vec.h"

using namespace std;



faust_core::faust_core() :
	data(std::vector<faust_spmat>()),
	factsMultiplied(false),
	factProduct(faust_mat()){}

faust_core(const std::vector<faust_spmat>& facts) :
	data(facts),
	factsMultiplied(false),
	factProduct(faust_mat()){}


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
	// factProduct = *(tmp[0]);
	// factsMultiplied = true;
}
#endif

