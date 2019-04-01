/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_Transform_GPU_HPP__
#define __FAUST_Transform_GPU_HPP__

#include "faust_Vect_gpu.h"
#include "faust_linear_algebra_gpu.h"
#include "faust_transform_algebra_gpu.h"
#include <iostream>
#include "faust_exception.h"
#include <fstream>
#include <iostream>
#include "cusparse.h"
#include "cublas_v2.h"

using namespace std;
template<typename FPP>
const char * Faust::Transform<FPP,Gpu>::m_className="Faust::Transform<FPP,Gpu>::";

template<typename FPP>
Faust::Transform<FPP,Gpu>::Transform() :
   data(std::vector<Faust::MatSparse<FPP,Gpu> >()),
   totalNonZeros(0)
{}


template<typename FPP>
Faust::Transform<FPP,Gpu>::Transform(const Faust::Transform<FPP,Gpu> & A) :
   data(A.data),
   totalNonZeros(A.totalNonZeros)
{}

template<typename FPP>
Faust::Transform<FPP,Gpu>::Transform(const std::vector<Faust::MatSparse<FPP,Gpu> >& facts, const FPP lambda_ /* =1.0 */) :
   data(facts),
   totalNonZeros(0)
{
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();

   if(lambda_ != 1.0 && data.size()>0)
      (data[0]) *= lambda_;



}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::get_facts(std::vector<Faust::MatDense<FPP,Gpu> >& facts)const
{
	facts.resize(size());
	for (int i=0;i<size();i++)
		facts[i] = data[i];

}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::print_file(const char* filename) const
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

template<typename FPP>
void Faust::Transform<FPP,Gpu>::scalarMultiply(const FPP scalar)
{
	if (size() > 0)
		data[0]*=scalar;
	else
		handleError(m_className,"scalarMultiply : empty faust can't be multiplied ");
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::init_from_file(const char* filename)
{
	clear();
	FILE* fp=fopen(filename,"r");
	int size_tmp;
	fscanf(fp,"%d\n", &size_tmp);
	fscanf(fp,"\n");
	Faust::MatSparse<FPP,Gpu> spmat;
	std::cout<<"size_tmp"<<size_tmp<<std::endl;
	for (int i=0;i<size_tmp;i++)
	{
		spmat.init_from_file(fp);
      Faust::MatSparse<FPP,Gpu> cu_spmat(spmat);
		data.push_back(cu_spmat);
	}
	fscanf(fp,"\n");
	updateNonZeros();
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::updateNonZeros()
{
	int totalNonZeros_tmp=0;
	for (int i=0;i<size();i++)
		totalNonZeros+=data[i].getNonZeros();
}


template<typename FPP>
Faust::Transform<FPP,Gpu>::Transform(const std::vector<Faust::MatDense<FPP,Gpu> >&facts) :
 data(std::vector<Faust::MatSparse<FPP,Gpu> >()),
   totalNonZeros(0)
{
	for (int i=0;i<facts.size();i++)
	{
		Faust::MatSparse<FPP,Gpu> spfact(facts[i]);
		data.push_back(spfact);
		totalNonZeros += data[i].getNonZeros();
	}
}



/*Faust::Transform<FPP,Gpu>::Transform(const Faust::Params& params) :
   data(std::vector<Faust::MatSparse<FPP,Gpu>>()),
   totalNonZeros(0)

{
   hierarchical_fact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   (data[0]) *= hier_fact.get_lambda();

}*/
template<typename FPP>
Faust::MatDense<FPP,Gpu> Faust::Transform<FPP,Gpu>::get_product(Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle)const
{
	//complexity of evaluating a Faust::Transform<FPP,Gpu>
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)
	if (size() == 0)
		handleError(m_className,"get_product : empty Faust::Transform<FPP,Gpu>");
	/*Faust::MatDense<FPP,Gpu> prod(data[0].getNbRow());

	if(getNbRow()<getNbCol())
	{
		prod.resize(getNbRow());
      cout<<"dim1="<<prod.getNbRow()<<" ; dim2="<<prod.getNbCol()<<endl;
		prod.setEyes();
		for(int i=0 ; i<data.size() ; i++)
		   //prod *= data[i];
         gemm(prod, data[i], prod, blasHandle, spblasHandle);
	}else
	{
		prod.resize(getNbCol());
      cout<<"dim1="<<prod.getNbRow()<<" ; dim2="<<prod.getNbCol()<<endl;
		prod.setEyes();
		for(int i=data.size()-1 ; i>=0 ; i--)
		   //prod.multiplyLeft(data[i]);
         gemm(data[i], prod, prod, spblasHandle);
	}*/

	Faust::MatDense<FPP,Gpu> prod;
	if ( (data[0].getNonZeros()+getNbRow()*(totalNonZeros-data[0].getNonZeros())) < (data[size()-1].getNonZeros()+getNbCol()*(totalNonZeros-data[size()-1].getNonZeros())) )
	{
		prod.init_from_cu_spmat(data[0], spblasHandle);
		for(int i=1 ; i<data.size() ; i++)
		{
		   //prod *= data[i];
         gemm(prod, data[i], prod, blasHandle, spblasHandle);
		}
	}else
	{
		prod.init_from_cu_spmat(data[size()-1], spblasHandle);
		for(int i=data.size()-2 ; i>=0 ; i--)
		{
		   //prod.multiplyLeft(data[i]);
         gemm(data[i], prod, prod, spblasHandle);
		}
	}



   return prod;
}

template<typename FPP>
int Faust::Transform<FPP,Gpu>::getNbRow() const
{
	//(size()!=0)?data[0].getNbRow():-1;
	if (size() !=0)
		return data[0].getNbRow();
	else
		return -1;

}
template<typename FPP>
int Faust::Transform<FPP,Gpu>::getNbCol() const
{
	//(size()!=0)?data[size()-1].getNbCol():-1;
	if (size() !=0)
		return data[size()-1].getNbCol();
	else
		return -1;
}

template<typename FPP>
FPP Faust::Transform<FPP,Gpu>::spectralNorm(const int nbr_iter_max, FPP threshold, int &flag) const
{
	if (size() == 0)
		return 1;
	else
	{
		Faust::Transform<FPP,Gpu> AtA(*this);
		AtA.transpose();
		if (getNbCol() < getNbRow())
			AtA.multiply(*this);
		else
			AtA.multiplyLeft(*this);
		return std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));
	}
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::multiply(const Faust::Transform<FPP,Gpu> & A)
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
				handleError(m_className,"multiply : dimensions of the 2 faustcore are in conflict");
			data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
		}
	}
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::multiplyLeft(const Faust::Transform<FPP,Gpu> & A)
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
				handleError(m_className,"multiplyLeft : dimensions of the 2 faustcore are in conflict");
			data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
		}
	}
}

template<typename FPP>
Faust::MatSparse<FPP,Gpu> Faust::Transform<FPP,Gpu>::get_fact(int id)const
{
	if((id>=size())||(id<0))
	{
		cout<<"id_fact error : "<<id<<endl;
		handleError(m_className,"get_fact : id exceed Faust::Transform<FPP,Gpu> size or id < 0");
	}
	cout<<"size_fact"<<size()<<endl;
	cout<<"id_fact"<<id<<endl;

	return data[id];
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::push_back(const Faust::MatSparse<FPP,Gpu>& S)
{
   if (size()>0)
      if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      handleError(m_className,"push_back : incorrect dimensions");
   data.push_back(S);
   totalNonZeros += S.getNonZeros();
}


template<typename FPP>
void Faust::Transform<FPP,Gpu>::push_first(const Faust::MatSparse<FPP,Gpu>& S)
{
   if (size()>0)
      if(data[0].getNbRow()!=S.getNbCol() || S.getNbRow()<1)
         handleError(m_className,"push_first : incorrect dimensions");
   data.insert(data.begin(),S);
   totalNonZeros += S.getNonZeros();
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::pop_back(Faust::MatSparse<FPP,Gpu>& S)
{
	if (size()>0)
	{
		S = data[size()-1];
		data.pop_back();
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("Faust::Transform<FPP,Gpu>::pop_back : empty Faust::Transform<FPP,Gpu>");
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::pop_first(Faust::MatSparse<FPP,Gpu>& S)
{
	if (size()>0)
	{
		S = data[0];
		data.erase(data.begin());
		totalNonZeros -= S.getNonZeros();
	}
	handleWarning("Faust::Transform<FPP,Gpu>::pop_back : empty Faust::Transform<FPP,Gpu>");
}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::pop_first(Faust::MatSparse<FPP,Gpu>& S) const
{
	if (size()>0)
	{
		S = data[0];
		//data.erase(data.begin());
		//totalNonZeros -= S.getNonZeros();
	}

}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::transpose()
{
	int nbFact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nbFact;i++)
		data[i].transpose();
}

//void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, Faust::SpBlasHandle<Gpu> spblasHandle);
template<typename FPP>
void Faust::Transform<FPP,Gpu>::mult( const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::SpBlasHandle<Gpu> spblasHandle)
{
	int nbFact=size();
	if (nbFact == 0)
	{
		cu_y=cu_x;
		handleWarning(m_className,"mult : empty Faust::Transform<FPP,Gpu>");

	}else
	{
		Faust::Vect<FPP,Gpu> vec_tmp(cu_x);
		for (int i=nbFact-1;i>=0;i--)
		{
			//std::cout<<"size x : "<< vec_tmp.size() <<" data nb col "<<data[i].getNbCol()<<std::endl;
			//std::cout<<"size y : "<< cu_y.size() <<" data nb row "<<data[i].getNbCol()<<std::endl;

			gemv(data[i],vec_tmp,cu_y,spblasHandle);
			vec_tmp=cu_y;
		}
	}

}

template<typename FPP>
void Faust::Transform<FPP,Gpu>::Display()const
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
