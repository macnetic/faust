#ifndef __FAUSTCORE_HPP__
#define __FAUSTCORE_HPP__

//#include "faust_core.h"
#include "faust_vec.h"
//#include "hierarchical_fact.h"
//#include "faust_params.h"
#include "LinAlgebra.h"
#include "faust_core_algebra.h"
#include <iostream>
#include "faust_exception.h"
#include <fstream>

using namespace std;
template<typename T>
const char * faust_core<T>::class_name="faust_core<T>::";

template<typename T>
faust_core<T>::faust_core() :
   data(std::vector<faust_spmat<T> >()),
   totalNonZeros(0)
{}

template<typename T>
faust_core<T>::faust_core(const faust_core<T> & A) :
   data(A.data),
   totalNonZeros(A.totalNonZeros)
{}

template<typename T>
faust_core<T>::faust_core(const std::vector<faust_spmat<T> > & facts, const T lambda_ /* =1.0 */) :
   data(facts),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();

   if(lambda_ != 1.0 && data.size()>0)
      (data[0]) *= lambda_;
}

template<typename T>
void faust_core<T>::get_facts(std::vector<faust_mat<T> >& facts)const
{
	facts.resize(size());
	for (int i=0;i<size();i++)
		facts[i] = data[i];
	
}




template<typename T>
void faust_core<T>::print_file(const char* filename) const
{
	if (size() > 0)
	{
		ofstream fichier;
		fichier.open(filename);	
		fichier<<size()<<endl<<endl;
		fichier.close();
		
		for (int i=0;i<size();i++)
		{
			data[i].print_file(filename,std::fstream::app);
		}
	}	
}

template<typename T>
void faust_core<T>::scalarMultiply(const T scalar)
{
	if (size() > 0)	
		data[0]*=scalar;
	else
		handleError(class_name,"scalarMultiply : empty faust can't be multiplied ");
}

template<typename T>
void faust_core<T>::init_from_file(const char* filename) 
{
	clear();
	FILE* fp=fopen(filename,"r");
	int size_tmp;
	fscanf(fp,"%d\n", &size_tmp);
	fscanf(fp,"\n");
	faust_spmat<T> spmat;
	std::cout<<"size_tmp"<<size_tmp<<std::endl;
	for (int i=0;i<size_tmp;i++)
	{
		spmat.init_from_file(fp);
		data.push_back(spmat);
	}
	fscanf(fp,"\n");
	updateNonZeros();
}

template<typename T>
void faust_core<T>::updateNonZeros()
{
	int totalNonZeros_tmp=0;
	for (int i=0;i<size();i++)
		totalNonZeros+=data[i].getNonZeros();
}


template<typename T>
faust_core<T>::faust_core(const std::vector<faust_mat<T> >&facts) :
 data(std::vector<faust_spmat<T> >()),
   totalNonZeros(0)
{	
	faust_spmat<T> spfact;
	for (int i=0;i<facts.size();i++)
	{
		spfact = facts[i];
		data.push_back(spfact);
		totalNonZeros += data[i].getNonZeros();
	}
}



/*faust_core<T>::faust_core(const faust_params& params) :
   data(std::vector<faust_spmat<T>>()),
   totalNonZeros(0)

{
   hierarchical_fact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   (data[0]) *= hier_fact.get_lambda();

}*/
template<typename T>
faust_mat<T> faust_core<T>::get_product()const
{
	//complexity of evaluating a faust_core 
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)	
	if (size() == 0)
	{
		handleError(class_name,"get_product : empty faust_core");						
	}	
	faust_mat<T> prod(data[0].getNbRow()); 
	
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
	/*faust_mat<T> prod;
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

template<typename T>
int faust_core<T>::getNbRow() const 
{
	if (size() != 0)
	{
		return data[0].getNbRow();
	}else
	{
		return -1;
	}
	
}
template<typename T>
int faust_core<T>::getNbCol() const 
{
	if (size() != 0)
	{
		return data[size()-1].getNbCol();
	}else
	{
		return -1;
	}
	
}

template<typename T>
T faust_core<T>::spectralNorm(const int nbr_iter_max, T threshold, int &flag) const
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

template<typename T>
void faust_core<T>::multiply(const faust_core & A)
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
				handleError(class_name,"multiply : dimensions of the 2 faustcore are in conflict");
			}
			data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
	}
	

template<typename T>	
void faust_core<T>::multiplyLeft(const faust_core & A)
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
				handleError(class_name,"multiplyLeft : dimensions of the 2 faustcore are in conflict");
			}
			data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
	}






template<typename T>
faust_spmat<T> faust_core<T>::get_fact(int id)const
{
	if((id>=size())||(id<0))
	{
		cout<<"id_fact error : "<<id<<endl;	
		handleError(class_name,"get_fact : id exceed faust_core size or id < 0"); 
	}
	cout<<"size_fact"<<size()<<endl;
	cout<<"id_fact"<<id<<endl;
	
	return data[id];
}

template<typename T>
void faust_core<T>::push_back(const faust_spmat<T>& S)
{
   if (size()>0)
   {   
      if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      {
		 handleError(class_name,"push_back : incorrect dimensions"); 
      }
   }
   data.push_back(S);
   totalNonZeros += S.getNonZeros();

	
}


template<typename T>
void faust_core<T>::push_first(const faust_spmat<T>& S)
{
   if (size()>0)
      if(data[0].getNbRow()!=S.getNbCol() || S.getNbRow()<1)
      {	
		handleError(class_name,"push_first : incorrect dimensions"); 
      }
   data.insert(data.begin(),S);
   totalNonZeros += S.getNonZeros();

	
}


template<typename T>
void faust_core<T>::pop_back(faust_spmat<T>& S)
{
	if (size()>0)
	{
		S = data[size()-1];
		data.pop_back();
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("faust_core<T>::pop_back : empty faust_core");	
}

template<typename T>
void faust_core<T>::pop_first(faust_spmat<T>& S)
{
	if (size()>0)
	{
		S = data[0];
		data.erase(data.begin());
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("faust_core<T>::pop_back : empty faust_core");		
}

template<typename T>
void faust_core<T>::pop_first(faust_spmat<T>& S) const 
{
	if (size()>0)
	{
		S = data[0];
		//data.erase(data.begin());
		//totalNonZeros -= S.getNonZeros();
	}
		
}

template<typename T>
void faust_core<T>::transpose()
{	
	int nb_fact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nb_fact;i++)
	{
		data[i].transpose();
	}
}






template<typename T>
void faust_core<T>::Display()const
{
   
	cout << "SIZE" << size() <<endl;
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

#endif
