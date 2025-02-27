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
#ifndef __FAUST_VECT_HPP__
#define __FAUST_VECT_HPP__

#include <iostream>
#include <iomanip>
#include <fstream>

template<typename FPP>
const char * Faust::Vect<FPP,Cpu>::m_className = "Faust::Vect<FPP,Cpu>::";

template<typename FPP>
template<typename FPP1>
Faust::Vect<FPP,Cpu>::Vect(const Faust::Vect<FPP1,Cpu>& v) : dim(v.dim), vec(v.dim)
{
for (int i=0;i<dim;i++)
	vec[i]=v(i);
}

template<typename FPP>
Faust::Vect<FPP,Cpu>::Vect(const Faust::Vect<FPP,Cpu>& v) : dim(v.dim), vec(v.vec)
{
}

template<typename FPP>
Faust::Vect<FPP,Cpu>::Vect(const faust_unsigned_int dim_, const FPP* data_) : dim(dim_), vec(dim_)
{
		memcpy(getData(), data_, dim*sizeof(FPP));
}

template<typename FPP>
Faust::Vect<FPP,Cpu>::Vect(Faust::Vect<FPP,Cpu>&& v)
{
	vec = std::move(v.vec);
	dim = v.dim;
}

template<typename FPP>
bool Faust::Vect<FPP,Cpu>::equality(Faust::Vect<FPP,Cpu> const &x, const Real<FPP> precision) const
{
  if (size() != x.size())
	handleError(m_className," equality : different dimension");

  for (int i=0;i<size();i++)
  {
	if (std::abs(x(i) - (*this)(i))>precision)
	{
		return false;
	}
  }

  return true;
}



template<typename FPP>
bool Faust::Vect<FPP,Cpu>::isReal() const
{
	bool isReal = (typeid(FPP) == typeid(double));
	isReal = (isReal || (typeid(FPP) == typeid(float)) );

	bool isComplex = (typeid(FPP) == typeid(std::complex<double>));
	isComplex = (isComplex || (typeid(FPP) == typeid(std::complex<float>)) );	
	
	if  ( (!isComplex) && (!isReal) )
	{
		handleError(m_className,"isReal : unknown type of scalar");
	}
	
		return isReal;
	
}






template<typename FPP>
void Faust::Vect<FPP,Cpu>::setOnes()
{
	vec.setOnes();
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::Display() const
{
    if(size()==0)
       std::cout << "empty vector";
    for (int i=0 ; i<size() ; i++)
	    std::cout << getData()[i] << " " ;
    std::cout<<std::endl<<std::endl;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::print_file(const char* filename)const
{
	std::ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<size() ; i++)
		fichier << std::setprecision(20) << vec (i) << std::endl;
	fichier.close();
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::resize(const int new_dim)
{

		if (new_dim <0)
		{
			handleError(m_className,"resize : new dimensions must be positive");
		}
		else if (dim != new_dim)
		{
			dim = new_dim;
			vec.conservativeResize(dim);

		}
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator*=(const FPP alpha)
{
   FPP*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] *= alpha;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator+=(const FPP alpha)
{
   FPP*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha;
}


template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator-=(const FPP alpha)
{
   FPP*const ptr_data = getData();
//#pragma omp parallel for
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] -= alpha;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator+=(const Faust::Vect<FPP,Cpu>& v)
{
   if(v.size()!=size())
   {

	  handleError(m_className,"operator+= : dimensions are in conflict");
   }
   vec += v.vec;
}


template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator-=(const Faust::Vect<FPP,Cpu>& v)
{
   if(v.size()!=size())
   {
	   handleError(m_className,"operator-= : dimensions are in conflict");
   }
	return *(this) -= v.getData();
}


template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator-=(const FPP* v_data)
{
   FPP*const ptr_data = getData();
//#pragma omp parallel for
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] -= v_data[i];
}

template<typename FPP>
FPP Faust::Vect<FPP,Cpu>::mean_relative_error(const Faust::Vect<FPP,Cpu>& v_ref)
{
   if(v_ref.size() != size())
   {
     handleError(m_className,"relative_error : sizes are different");

   }

   Faust::Vect<FPP,Cpu> tmp(size());
   for(int i=0 ; i<size() ; i++)
	  if(vec[i] != FPP(0) && v_ref.vec[i] != FPP(0))
		  tmp[i] = Faust::fabs((vec[i]-v_ref.vec[i])/v_ref.vec[i]);
	  else
		  tmp[i] = 0;

   return tmp.mean();
}

template<typename FPP>
FPP Faust::Vect<FPP,Cpu>::mean()
{
	FPP m = 0;
	for(int i=0; i < size(); i++)
		m += getData()[i];
	m /= size();
	return m;
}


template<typename FPP>
template<typename FPP1>
void Faust::Vect<FPP,Cpu>::operator=(Faust::Vect<FPP1,Cpu> const& y)
{

	resize(y.dim);
//#pragma omp parallel for
	for (int i=0;i<dim;i++)
		vec[i]=y(i);

}

template<typename FPP>
Faust::Vect<FPP,Cpu>& Faust::Vect<FPP,Cpu>::operator=(const Faust::Vect<FPP,Cpu>& y)
{
	if(this->size() != y.size())
		this->resize(y.size());
	memcpy(this->getData(), y.getData(), y.size()*sizeof(FPP));
	return *this;
}

template<typename FPP>
Faust::Vect<FPP,Cpu>& Faust::Vect<FPP,Cpu>::operator=(Faust::Vect<FPP,Cpu>&& y)
{
	this->vec = std::move(y.vec);
	this->dim = y.dim;
	return *this;
}

template <typename FPP>
FPP Faust::Vect<FPP,Cpu>::dot(const Faust::Vect<FPP,Cpu>& v) const
{
//return dot(*this, v);
   if(size() != v.size())
      handleError("linear_algebra","dot : the two vectors don't have the same size");
   FPP result = vec.dot(v.vec);
   return result;

}



template<typename FPP>
void  Faust::Vect<FPP,Cpu>::multiplyLeft(Faust::MatDense<FPP,Cpu> const& A)
{
   Faust::gemv<FPP>(A, *this, *this, 1.0, 0.0, 'N');
}

template<typename FPP>
void  Faust::Vect<FPP,Cpu>::multiplyLeft(Faust::MatSparse<FPP,Cpu> const& S,const char transS)
{
	faust_unsigned_int nbColOpS,nbRowOpS;
	S.setOp(transS,nbRowOpS,nbColOpS);

	if(nbColOpS != vec.size())
	{
		 handleError(m_className,"multiplyLeft : incorrect dimensions");

	}

		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
		#endif
		if (transS == 'N')
			vec = S.mat * vec;
		else if(transS == 'T')
		        vec = S.mat.transpose() * vec;
		else if(transS == 'H')
		        vec = S.mat.adjoint() * vec;
		else
			throw std::runtime_error("Unknown op transS");

		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			//std::cout <<"0 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<std::endl;
			t_local_multiplyLeft.reset();
		#endif


	this->dim = nbRowOpS;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::conjugate()
{
	vec = vec.conjugate().eval();
}

template<typename FPP>
FPP Faust::Vect<FPP, Cpu>::normL1() const
{
	faust_unsigned_int j;
	FPP sum = 0;
	for(j=0;j<this->size();j++)
		sum += std::abs(vec.data()[j]);
	return sum;
}

	template<typename FPP>
Faust::Vect<FPP, Cpu>* Faust::Vect<FPP, Cpu>::rand(faust_unsigned_int size)
{
	Faust::Vect<FPP, Cpu>* v = nullptr;
	try
	{
		v = new Faust::Vect<FPP, Cpu>(size);
		v->setRand();
	}
	catch(std::bad_alloc e)
	{
		//		std::cerr << "Out of memory." << std::endl;
	}
	return v;
}

	template<typename FPP>
void Faust::Vect<FPP, Cpu>::setRand()
{
	vec = Eigen::Matrix<FPP, Eigen::Dynamic, 1>::Random(this->size(), 1);
}

template<typename FPP>
FPP Faust::Vect<FPP, Cpu>::max_coeff(int *index) const
{
	FPP max = getData()[0]; //FPP(std::numeric_limits<double>::max());
	FPP e;
	*index = 0;
	//			vec.getData()[i] = Eigen::abs(mat.row(i)).maxCoeff(col_indices+i);
	for(int j=0;j<this->size(); j++)
	{
		e = getData()[j];
		if(Faust::fabs(e) > Faust::fabs(max))
		{
			max = e;
			*index = j;
		}
	}
	return max;
}

template<typename FPP>
void Faust::Vect<FPP, Cpu>::normalize()
{
	auto n = norm();
	if(n == Real<FPP>(0))
		throw std::domain_error("Can't normalize a zero-norm vector.");
	scalarMultiply(1/norm());
}

#endif
