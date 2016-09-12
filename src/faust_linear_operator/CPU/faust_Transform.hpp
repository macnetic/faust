/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __FAUSTCORE_HPP__
#define __FAUSTCORE_HPP__


#include "faust_Vect.h"
//#include "faust_HierarchicalFact.h"
//#include "faust_Params.h"
#include "faust_linear_algebra.h"
#include "faust_transform_algebra.h"
#include <iostream>
#include "faust_exception.h"
#include <fstream>
#include "faust_BlasHandle.h"
#include "faust_SpBlasHandle.h"





template<typename FPP>
void Faust::Transform<FPP,Cpu>::faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB) const
{

	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if (size() == 0)
	{
		handleError(m_className,"faust_gemm : empty Faust::Transform");
	}
	if (typeA == 'T')
	{
		nbRowOpA = this->getNbCol();
		nbColOpA = this->getNbRow();
	}else
	{
		nbRowOpA = this->getNbRow();
		nbColOpA = this->getNbCol();
	}


	if (typeB == 'T')
	{
		nbRowOpB = B.getNbCol();
		nbColOpB = B.getNbRow();
	}else
	{
		nbRowOpB = B.getNbRow();
		nbColOpB = B.getNbCol();
	}


	if (nbColOpA != nbRowOpB)
	{
		handleError(m_className, "faust_gemm : dimension conflict  between matrix op((*this)) and matrix op(B)");

	}

    int* ind_ptr = new int[size()];
    for (int j=0 ; j<size(); j++)
    {
    	if (typeA == 'T')
        	ind_ptr[j] = j;
    	else
        	ind_ptr[j] =size()-1-j;
    }

    if (beta == 0)
    {
	data[ind_ptr[0]].faust_gemm(B,C,alpha,beta,typeA,typeB);
	if (size() > 1)
	{
		Faust::MatDense<FPP,Cpu> tmp1(C);
    		for (int i=1;i<size();i++)
    		{
			data[ind_ptr[i]].faust_gemm(tmp1,C,1.0,0.0,typeA,'N');
			tmp1=C;
		}
    	}

    }else
    {
	if (( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		handleError(m_className, "faust_gemm : invalid dimension for output matrix C");
	}


	if (size() ==1)
	{
		data[0].faust_gemm(B,C,alpha,beta,typeA,typeB);
	}else
	{
		Faust::MatDense<FPP,Cpu> tmp2,tmp1;
		data[ind_ptr[0]].faust_gemm(B,tmp1,alpha,0.0,typeA,typeB);

		for (int i=1;i<(size()-1);i++)
    		{
			data[ind_ptr[i]].faust_gemm(tmp1,tmp2,1.0,0.0,typeA,'N');
			tmp1=tmp2;
		}

    		data[ind_ptr[size()-1]].faust_gemm(tmp1,C,1.0,beta,typeA,'N');
	}

   }
   delete[] ind_ptr;
    ind_ptr = NULL;

}


template<typename FPP>
const char * Faust::Transform<FPP,Cpu>::m_className="Faust::Transform<FPP,Cpu>";


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform() :
	data(std::vector<Faust::MatSparse<FPP,Cpu> >()),
	totalNonZeros(0)
{

#ifdef __COMPILE_TIMERS__
	t_multiply_vector.resize(0);
#endif 

}


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const Faust::Transform<FPP,Cpu> & A) :
	data(A.data),
	totalNonZeros(A.totalNonZeros)
{
#ifdef __COMPILE_TIMERS__
		this->t_multiply_vector.resize(data.size());
#endif 

}

template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const std::vector<Faust::MatSparse<FPP,Cpu> > & facts, const FPP lambda_ /* =1.0 */) :
    data(facts),
	totalNonZeros(0)
{
	for (int i=0 ; i<data.size() ; i++)
		totalNonZeros += data[i].getNonZeros();

	if(lambda_ != 1.0 && data.size()>0)
		(data[0]) *= lambda_;
	
	#ifdef __COMPILE_TIMERS__
		this->t_multiply_vector.resize(data.size());
	#endif
	this->check_factors_validity();
}

template<typename FPP>

void Faust::Transform<FPP,Cpu>::check_factors_validity() const
{
		
	if (size() > 0)
	{
	 
	 for (int i=0;i<=size()-2;i++)
	 {
	    if (data[i].getNbCol() != data[i+1].getNbRow())
		handleError(m_className,"check_factors_validity : dimensions of the factors mismatch");
			
	 }	 	


	}

}


template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_facts(std::vector<Faust::MatDense<FPP,Cpu> >& facts)const
{
	facts.resize(size());
	for (int i=0;i<size();i++)
		facts[i] = data[i];
}



template<typename FPP>
void Faust::Transform<FPP,Cpu>::print_file(const char* filename) const
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

template<typename FPP>
void Faust::Transform<FPP,Cpu>::scalarMultiply(const FPP scalar)
{
	if (size() > 0)
		data[0]*=scalar;
	else
		handleError(m_className,"scalarMultiply : empty faust can't be multiplied ");
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::init_from_file(const char* filename)
{
	clear();
	FILE* fp=fopen(filename,"r");
	int size_tmp;
	fscanf(fp,"%d\n", &size_tmp);
	fscanf(fp,"\n");
	Faust::MatSparse<FPP,Cpu> spmat;
	std::cout<<"size_tmp"<<size_tmp<<std::endl;
	for (int i=0;i<size_tmp;i++)
	{
		spmat.init_from_file(fp);
		data.push_back(spmat);
	}
	fscanf(fp,"\n");
	updateNonZeros();
	
	#ifdef __COMPILE_TIMERS__
		this.t_multiply_vector.resize(data.size());
	#endif
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::updateNonZeros()
{
	int totalNonZeros_tmp=0;
	for (int i=0;i<size();i++)
		totalNonZeros+=data[i].getNonZeros();
}


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const std::vector<Faust::MatDense<FPP,Cpu> >&facts) :
	data(std::vector<Faust::MatSparse<FPP,Cpu> >()),
	totalNonZeros(0)
{
	Faust::MatSparse<FPP,Cpu> spfact;
	for (int i=0;i<facts.size();i++)
	{
		spfact = facts[i];
		data.push_back(spfact);
		totalNonZeros += data[i].getNonZeros();
	}
	#ifdef __COMPILE_TIMERS__
		this.t_multiply_vector.resize(data.size());
	#endif
}



/*Faust::Transform<FPP,Cpu>::Transform(const Faust::Params& params) :
   data(std::vector<Faust::MatSparse<FPP,Cpu>>()),
   totalNonZeros(0)

{
   Faust::HierarchicalFact hier_fact(params);
   hier_fact.compute_facts();
   hier_fact.get_facts(data);
   for (int i=0 ; i<data.size() ; i++)
      totalNonZeros += data[i].getNonZeros();
   (data[0]) *= hier_fact.get_lambda();

}*/
template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::get_product(const char opThis)const
{
	//complexity of evaluating a Faust::Transform
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)

	if (size() == 0)
	{
		handleError(m_className,"get_product : empty Faust::Transform");
	}
	
	Faust::MatDense<FPP,Cpu> prod(data[0].getNbRow());
	
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
	if (opThis == 'T')
		prod.transpose();

	/*Faust::MatDense<FPP,Cpu> prod;
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

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::get_product(Faust::BlasHandle<Cpu> blas_handle, Faust::SpBlasHandle<Cpu> spblas_handle)const
{
    return (*this).get_product();
}


template<typename FPP>
faust_unsigned_int Faust::Transform<FPP,Cpu>::getNbRow() const
{
	if (size() != 0)
	{
		return data[0].getNbRow();
	}else
	{
		return 0;
	}

}

template<typename FPP>
faust_unsigned_int Faust::Transform<FPP,Cpu>::getNbCol() const
{
	if (size() != 0)
	{
		return data[size()-1].getNbCol();
	}else
	{
		return 0;
	}
}

template<typename FPP>
FPP Faust::Transform<FPP,Cpu>::spectralNorm(const int nbr_iter_max, FPP threshold, int &flag) const
{
	if (size() == 0)
	{
		return 1;
	}else
	{
		Faust::Transform<FPP,Cpu> AtA((*this)); // modif AL <FPP,Cpu>
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

template<typename FPP>
void Faust::Transform<FPP,Cpu>::multiply(const Faust::Transform<FPP,Cpu> & A) // modif AL <FPP,Cpu>
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
				handleError(m_className,"multiply : dimensions of the 2 faust_transform are in conflict");
			}
			data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
		}
	}
	#ifdef __COMPILE_TIMERS__
		this.t_multiply_vector.resize(data.size());
	#endif
}


template<typename FPP>
void Faust::Transform<FPP,Cpu>::multiplyLeft(const Faust::Transform<FPP,Cpu> & A) // modif AL <FPP,Cpu>
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
				handleError(m_className,"multiplyLeft : dimensions of the 2 faustcore are in conflict");
			}
			data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
		}
	}
}



template<typename FPP>
Faust::MatSparse<FPP,Cpu> Faust::Transform<FPP,Cpu>::get_fact(faust_unsigned_int id)const
{
	if(id>=size())
	{
		handleError(m_className,"get_fact : id exceed Faust::Transform size or id < 0");
	}


	return data[id];
}







template<typename FPP>
void Faust::Transform<FPP,Cpu>::push_back(const Faust::MatSparse<FPP,Cpu>& S)
{
	if (size()>0)
	{
		if(data[size()-1].getNbCol()!=S.getNbRow() || S.getNbRow()<1)
      		{
			handleError(m_className,"push_back : incorrect dimensions");
     		}
   	}
   	data.push_back(S);
   	totalNonZeros += S.getNonZeros();
	#ifdef __COMPILE_TIMERS__
		this->t_multiply_vector.push_back(Faust::Timer());
	#endif
	
}


template<typename FPP>
void Faust::Transform<FPP,Cpu>::push_first(const Faust::MatSparse<FPP,Cpu>& S)
{
	if (size()>0)
		if(data[0].getNbRow()!=S.getNbCol() || S.getNbRow()<1)
      		{
			handleError(m_className,"push_first : incorrect dimensions");
      		}
	data.insert(data.begin(),S);
   	totalNonZeros += S.getNonZeros();
	
	#ifdef __COMPILE_TIMERS__
		this->t_multiply_vector.insert(this->t_multiply_vector().begin(),Faust::Timer());
	#endif
}


template<typename FPP>
void Faust::Transform<FPP,Cpu>::pop_back(Faust::MatSparse<FPP,Cpu>& S)
{
	if (size()>0)
	{
		S = data[size()-1];
		data.pop_back();
		totalNonZeros -= S.getNonZeros();
		#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector.pop_back();
		#endif
	}
	handleWarning("Faust::Transform<FPP,Cpu>::pop_back : empty Faust::Transform");
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::pop_first(Faust::MatSparse<FPP,Cpu>& S)
{
	if (size()>0)
	{
		S = data[0];
		data.erase(data.begin());
		totalNonZeros -= S.getNonZeros();
		#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector.erase(this->t_multiply_vector.begin());
		#endif
	}
	handleWarning("Faust::Transform<FPP,Cpu>::pop_back : empty Faust::Transform");
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::pop_first(Faust::MatSparse<FPP,Cpu>& S) const
{
	if (size()>0)
	{
		S = data[0];
		//data.erase(data.begin());
		//totalNonZeros -= S.getNonZeros();
	}

}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::transpose()
{
	int nbFact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nbFact;i++)
	{
		data[i].transpose();
	}
}
template<typename FPP>
#ifdef __COMPILE_TIMERS__
	Faust::Vect<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::Vect<FPP,Cpu> x,const char opThis)
#else
	Faust::Vect<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::Vect<FPP,Cpu> x,const char opThis) const
#endif
{
	#ifdef __COMPILE_TIMERS__
			if (this->t_multiply_vector.size() != this->data.size())
				handleError(m_className,"multiply : invalid t_multiply_vector size to measure time");
		#endif

	if (size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> : multiply : empty Faust::Transform<FPP,Cpu>");
	
	Faust::Vect<FPP,Cpu> vec(x);
 
	if (opThis == 'N')
	{	
		for (int i=this->size()-1 ; i >= 0 ; i--)
		{
			#ifdef __COMPILE_TIMERS__
				this->t_multiply_vector[i].start();
			#endif

			vec.multiplyLeft(data[i]);

			#ifdef __COMPILE_TIMERS__
				this->t_multiply_vector[i].stop();
			#endif
		}
	}else
	{		
		for (int i=0 ; i < this->size() ; i++)
		{	
			#ifdef __COMPILE_TIMERS__
				this->t_multiply_vector[i].start();
			#endif

			vec.multiplyLeft(data[i],opThis);

			#ifdef __COMPILE_TIMERS__
				this->t_multiply_vector[i].stop();
			#endif
		}
		
	}
	return vec;

		
	
}




template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::MatDense<FPP,Cpu> A,const char opThis) const
{


	if (size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> : multiply : empty Faust::Transform<FPP,Cpu>");
	
	Faust::MatDense<FPP,Cpu> mat(A);
 
	if (opThis == 'N')
	{	
		for (int i=this->size()-1 ; i >= 0 ; i--)
			mat.multiplyLeft(data[i]);
		
	}else
	{		
		for (int i=0 ; i < this->size() ; i++)
		{	
			mat.multiplyLeft(data[i],opThis);
		}
		
	}
	return mat;

		
	
}








template<typename FPP>
void Faust::Transform<FPP,Cpu>::setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const
{
	if (size() > 0)
	{	
		if(op == 'N')
    		{
        		nbRowOp=data[0].getNbRow();
        		nbColOp=data[size()-1].getNbCol();
    		}
    		else if(op == 'T')
    		{
        		nbRowOp=data[size()-1].getNbCol();
        		nbColOp=data[0].getNbRow();
    		}
    		else
        		handleError(m_className,"setOp : invalid character");
	}else
	{
	   nbRowOp = 0;
	   nbColOp = 0;	
	   handleWarning("Faust::Transform<FPP,Cpu>::setOp : empty Faust::Transform");
	}
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::Display()const
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







#ifdef __COMPILE_TIMERS__
	template<typename FPP>
	void Faust::Transform<FPP,Cpu>::print_timers() const
	{
		std::cout << "timers in Faust::Transform :" << endl;
		std::cout << "t_multiply_vector :" << endl;	
		float sum_tps=0;
		for (int i=0;i<t_multiply_vector.size();i++)
			sum_tps+=t_multiply_vector[i].get_time();		
		float current_time;		
		for (int i=0;i<t_multiply_vector.size();i++)
		{
			float current_time=t_multiply_vector[i].get_time();
			std::cout<<"	FACTOR "<<i<<", size "<<data[i].getNbRow()<<"x"<<data[i].getNbCol()<<" ,nnz "<<data[i].getNonZeros()<<std::endl;
			std::cout<<" 	,tps "<<current_time<<" ,% "<<100*current_time/sum_tps<< " ,nb_calls  "<<t_multiply_vector[i].get_nb_call()<<std::endl;
		}
		std::cout<<"	total tps "<<sum_tps<<std::endl; 
	


	}
#endif








#endif

