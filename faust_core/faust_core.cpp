#include "faust_core.h"
#include "faust_vec.h"
#include "hierarchical_fact.h"
#include "faust_params.h"
#include "LinAlgebra.h"

using namespace std;



faust_core::faust_core() :
   data(std::vector<faust_spmat>()),
   totalNonZeros(0)
{}


faust_core::faust_core(const faust_core & A) :
   data(A.data),
   totalNonZeros(A.totalNonZeros)
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

faust_mat faust_core::get_product()const
{
	//complexity of evaluating a faust_core 
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)	
	if (size() == 0)
	{
		cerr << "Error in faust_core::get_product : empty_faust_core" << endl;
				exit(EXIT_FAILURE);
	}	
	faust_mat prod(data[0].getNbRow()); 
	
	if(getNbRow()<getNbCol())
	{
		prod.resize(getNbRow());
		prod.setEyes();
		for(int i=0 ; i<data.size() ; i++)
		prod *= data[i];		
	}else
	{
		prod.resize(getNbCol());
		prod.setEyes();
		for(int i=data.size()-1 ; i>=0 ; i--)
		prod.multiplyLeft(data[i]);		
	}
	/*faust_mat prod;
	if ( (data[0].getNonZeros()+getNbRow()*(totalNonZeros-data[0].getNonZeros())) < (data[size()-1].getNonZeros()+getNbCol()*(totalNonZeros-data[size()-1].getNonZeros())) )
	{
		prod = data[0];	
		for(int i=1 ; i<data.size() ; i++)
		prod *= data[i];	
	}else
	{	
		prod = data[size()-1];	
		for(int i=data.size()-2 ; i>=0 ; i--)
		prod.multiplyLeft(data[i]);	
	}	*/
   
   

   return prod;
}


int faust_core::getNbRow() const 
{
	if (size() != 0)
	{
		return data[0].getNbRow();
	}else
	{
		return -1;
	}
	
}

int faust_core::getNbCol() const 
{
	if (size() != 0)
	{
		return data[size()-1].getNbCol();
	}else
	{
		return -1;
	}
	
}





faust_real faust_core::spectralNorm(const int nbr_iter_max, faust_real threshold, int &flag) const
{
	if (size() == 0)
	{
		return 1;
	}else
	{
		faust_core AtA((*this));
		AtA.transpose();
		if (getNbCol() < getNbRow())
		{
			
			AtA.multiply((*this));
			
		}else
		{
			AtA.multiplyLeft((*this));	
		}		
		return std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));
		
		
		
	}	
	
	
	
}












void faust_core::multiply(const faust_core & A)
{
	if (A.size() == 0)
	{
		
	}else
	{
		if (size() == 0)
		{
			(*this)=A;
		}
		else
		{
			if (getNbCol() != A.getNbRow())
			{
				cerr << "Error in faust_core::multiply : dimensions of the 2 faustcore are in conflict" << endl;
				exit(EXIT_FAILURE);
			}
			data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
	}
	
	
void faust_core::multiplyLeft(const faust_core & A)
{
	if (A.size() == 0)
	{
		
	}else
	{
		if (size() == 0)
		{
			(*this)=A;
		}
		else
		{
			if (getNbRow() != A.getNbCol())
			{
				cerr << "Error in faust_core::multiplyLeft : dimensions of the 2 faustcore are in conflict" << endl;
				exit(EXIT_FAILURE);
			}
			data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
	}







faust_spmat faust_core::get_fact(int id)const
{
	if(id>=size())
	{
		cerr << "Error in faust_core::get_fact : id exceed faust_core size" << endl;
         exit(EXIT_FAILURE);
	}
	return data[id];
}


void faust_core::push_back(const faust_spmat& S)
{
   if (size()>0)
   {   
      if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      {
         cerr << "Error in faust_core::push_back : incorrect dimensions" << endl;
         exit(EXIT_FAILURE);
      }
   }
   data.push_back(S);
   totalNonZeros += S.getNonZeros();

	
}



void faust_core::push_first(const faust_spmat& S)
{
   if (size()>0)
      if(data[0].getNbRow()!=S.getNbCol() || S.getNbRow()<1)
      {
         cerr << "Error in faust_core::push_first : incorrect dimensions" << endl;
         exit(EXIT_FAILURE);
      }
   data.emplace(data.begin(),S);
   totalNonZeros += S.getNonZeros();

	
}




void faust_core::pop_back(faust_spmat& S)
{
	if (size()>0)
	{
		S = data[size()-1];
		data.pop_back();
		totalNonZeros -= S.getNonZeros();
	}
		
}


void faust_core::pop_first(faust_spmat& S)
{
	if (size()>0)
	{
		S = data[0];
		data.erase(data.begin());
		totalNonZeros -= S.getNonZeros();
	}
		
}


void faust_core::pop_first(faust_spmat& S) const 
{
	if (size()>0)
	{
		S = data[0];
		//data.erase(data.begin());
		//totalNonZeros -= S.getNonZeros();
	}
		
}





void faust_core::transpose()
{	
	int nb_fact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nb_fact;i++)
	{
		data[i].transpose();
	}
}




/*void faust_core::operator*=(const faust_core&  f)
{
   for (int i=0 ; i<f.size() ; i++)
      push_back(f.data[i]);
}*/


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
