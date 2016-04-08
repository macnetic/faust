#ifndef __FAUST_CORE_CU_HPP__
#define __FAUST_CORE_CU_HPP__
//#include "faust_core_cu.h"
#include "faust_cu_vec.h"
//#include "hierarchical_fact.h"
//#include "faust_params.h"
#include "LinAlgebra_cu.h"
#include "faust_core_algebra_cu.h"
#include <iostream>
#include "faust_exception.h"
#include <fstream>
#include <iostream>
using namespace std;
template<typename T>
const char * faust_core_cu<T>::class_name="faust_core_cu<T>::";

template<typename T>
faust_core_cu<T>::faust_core_cu() :
   data(std::vector<faust_cu_spmat<T> >()),
   totalNonZeros(0)
{}


template<typename T>
faust_core_cu<T>::faust_core_cu(const faust_core_cu<T> & A) :
   data(A.data),
   totalNonZeros(A.totalNonZeros)
{}

template<typename T>
faust_core_cu<T>::faust_core_cu(const std::vector<faust_cu_spmat<T> >& facts, const T lambda_ /* =1.0 */) :
   data(facts),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();

   if(lambda_ != 1.0 && data.size()>0)
      (data[0]) *= lambda_;

  
   
}

template<typename T>
void faust_core_cu<T>::get_facts(std::vector<faust_cu_mat<T> >& facts)const
{
	facts.resize(size());
	for (int i=0;i<size();i++)
		facts[i] = data[i];
	
}

template<typename T>
void faust_core_cu<T>::print_file(const char* filename) const
{
	if (size() > 0)
	{
		ofstream fichier;
		fichier.open(filename);	
		fichier<<size()<<endl<<endl;
		fichier.close();
		
		for (int i=0;i<size();i++)
			data[i].print_file(filename,std::fstream::app);
	}	
}

template<typename T>
void faust_core_cu<T>::scalarMultiply(const T scalar)
{
	if (size() > 0)	
		data[0]*=scalar;
	else
		handleError(class_name,"scalarMultiply : empty faust can't be multiplied ");
}

template<typename T>
void faust_core_cu<T>::init_from_file(const char* filename) 
{
	clear();
	FILE* fp=fopen(filename,"r");
	int size_tmp;
	fscanf(fp,"%d\n", &size_tmp);
	fscanf(fp,"\n");
	faust_cu_spmat<T> spmat;
	std::cout<<"size_tmp"<<size_tmp<<std::endl;
	for (int i=0;i<size_tmp;i++)
	{
		spmat.init_from_file(fp);
      faust_cu_spmat<T> cu_spmat(spmat);
		data.push_back(cu_spmat);
	}
	fscanf(fp,"\n");
	updateNonZeros();
}

template<typename T>
void faust_core_cu<T>::updateNonZeros()
{
	int totalNonZeros_tmp=0;
	for (int i=0;i<size();i++)
		totalNonZeros+=data[i].getNonZeros();
}


template<typename T>
faust_core_cu<T>::faust_core_cu(const std::vector<faust_cu_mat<T> >&facts) :
 data(std::vector<faust_cu_spmat<T> >()),
   totalNonZeros(0)
{	
	for (int i=0;i<facts.size();i++)
	{
		faust_cu_spmat<T> spfact(facts[i]);
		data.push_back(spfact);
		totalNonZeros += data[i].getNonZeros();
	}
}



/*faust_core_cu<T>::faust_core_cu(const faust_params& params) :
   data(std::vector<faust_cu_spmat<T>>()),
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
faust_cu_mat<T> faust_core_cu<T>::get_product(cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)const
{
	//complexity of evaluating a faust_core_cu 
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)	
	if (size() == 0)
		handleError(class_name,"get_product : empty faust_core_cu");						
	/*faust_cu_mat<T> prod(data[0].getNbRow());

	if(getNbRow()<getNbCol())
	{
		prod.resize(getNbRow());
      cout<<"dim1="<<prod.getNbRow()<<" ; dim2="<<prod.getNbCol()<<endl;
		prod.setEyes();
		for(int i=0 ; i<data.size() ; i++)
		   //prod *= data[i];		
         gemm(prod, data[i], prod, cublasHandle, cusparseHandle);
	}else
	{
		prod.resize(getNbCol());
      cout<<"dim1="<<prod.getNbRow()<<" ; dim2="<<prod.getNbCol()<<endl;
		prod.setEyes();
		for(int i=data.size()-1 ; i>=0 ; i--)
		   //prod.multiplyLeft(data[i]);		
         gemm(data[i], prod, prod, cusparseHandle);
	}*/

	faust_cu_mat<T> prod;
	if ( (data[0].getNonZeros()+getNbRow()*(totalNonZeros-data[0].getNonZeros())) < (data[size()-1].getNonZeros()+getNbCol()*(totalNonZeros-data[size()-1].getNonZeros())) )
	{
		prod.init_from_cu_spmat(data[0], cusparseHandle);	
		for(int i=1 ; i<data.size() ; i++)
		{
		   //prod *= data[i];	
         gemm(prod, data[i], prod, cublasHandle, cusparseHandle);
		}
	}else
	{	
		prod.init_from_cu_spmat(data[size()-1], cusparseHandle);	
		for(int i=data.size()-2 ; i>=0 ; i--)
		{
		   //prod.multiplyLeft(data[i]);	
         gemm(data[i], prod, prod, cusparseHandle);
		}
	}	
   
   

   return prod;
}

template<typename T>
int faust_core_cu<T>::getNbRow() const 
{
	//(size()!=0)?data[0].getNbRow():-1;
	if (size() !=0)
		return data[0].getNbRow();
	else
		return -1;
	
}
template<typename T>
int faust_core_cu<T>::getNbCol() const 
{
	//(size()!=0)?data[size()-1].getNbCol():-1;
	if (size() !=0)
		return data[size()-1].getNbCol();
	else
		return -1;	
}

template<typename T>
T faust_core_cu<T>::spectralNorm(const int nbr_iter_max, T threshold, int &flag) const
{
	if (size() == 0)
		return 1;
	else
	{
		faust_core_cu AtA(*this);
		AtA.transpose();
		if (getNbCol() < getNbRow())
			AtA.multiply(*this);
		else
			AtA.multiplyLeft(*this);	
		return std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));	
	}	
}

template<typename T>
void faust_core_cu<T>::multiply(const faust_core_cu & A)
{
	if (A.size() == 0)
	{	
	}
   else
	{
		if (size() == 0)
			(*this)=A;
		else
		{
			if (getNbCol() != A.getNbRow())
				handleError(class_name,"multiply : dimensions of the 2 faustcore are in conflict");
			data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
}
	
template<typename T>	
void faust_core_cu<T>::multiplyLeft(const faust_core_cu & A)
{
	if (A.size() == 0)
	{		
	}
   else
	{
		if (size() == 0)
			(*this)=A;
		else
		{
			if (getNbRow() != A.getNbCol())
				handleError(class_name,"multiplyLeft : dimensions of the 2 faustcore are in conflict");
			data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;	
		}
	}	
}

template<typename T>
faust_cu_spmat<T> faust_core_cu<T>::get_fact(int id)const
{
	if((id>=size())||(id<0))
	{
		cout<<"id_fact error : "<<id<<endl;	
		handleError(class_name,"get_fact : id exceed faust_core_cu size or id < 0"); 
	}
	cout<<"size_fact"<<size()<<endl;
	cout<<"id_fact"<<id<<endl;
	
	return data[id];
}

template<typename T>
void faust_core_cu<T>::push_back(const faust_cu_spmat<T>& S)
{
   if (size()>0)
      if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      handleError(class_name,"push_back : incorrect dimensions"); 
   data.push_back(S);
   totalNonZeros += S.getNonZeros();
}


template<typename T>
void faust_core_cu<T>::push_first(const faust_cu_spmat<T>& S)
{
   if (size()>0)
      if(data[0].getNbRow()!=S.getNbCol() || S.getNbRow()<1)
         handleError(class_name,"push_first : incorrect dimensions"); 
   data.insert(data.begin(),S);
   totalNonZeros += S.getNonZeros();
}

template<typename T>
void faust_core_cu<T>::pop_back(faust_cu_spmat<T>& S)
{
	if (size()>0)
	{
		S = data[size()-1];
		data.pop_back();
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("faust_core_cu<T>::pop_back : empty faust_core_cu");	
}

template<typename T>
void faust_core_cu<T>::pop_first(faust_cu_spmat<T>& S)
{
	if (size()>0)
	{
		S = data[0];
		data.erase(data.begin());
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("faust_core_cu<T>::pop_back : empty faust_core_cu");		
}

template<typename T>
void faust_core_cu<T>::pop_first(faust_cu_spmat<T>& S) const 
{
	if (size()>0)
	{
		S = data[0];
		//data.erase(data.begin());
		//totalNonZeros -= S.getNonZeros();
	}
		
}

template<typename T>
void faust_core_cu<T>::transpose()
{	
	int nb_fact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nb_fact;i++)
		data[i].transpose();
}

//void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cusparseHandle_t cusparseHandle);
template<typename T>
void faust_core_cu<T>::mult( const faust_cu_vec<T>& cu_x, faust_cu_vec<T>& cu_y, cusparseHandle_t cusparseHandle)
{
	int nb_fact=size();
	if (nb_fact == 0)
	{	
		cu_y=cu_x;
		handleWarning(class_name,"mult : empty faust_core_cu");

	}else
	{
		faust_cu_vec<T> vec_tmp(cu_x);
		for (int i=nb_fact-1;i>=0;i--)
		{	
			//std::cout<<"size x : "<< vec_tmp.size() <<" data nb col "<<data[i].getNbCol()<<std::endl;
			//std::cout<<"size y : "<< cu_y.size() <<" data nb row "<<data[i].getNbCol()<<std::endl;

			gemv(data[i],vec_tmp,cu_y,cusparseHandle);
			vec_tmp=cu_y;	
		}
	}

}

template<typename T>
void faust_core_cu<T>::Display()const
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
