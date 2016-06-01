#ifndef __FAUST_VECT_HPP__
#define __FAUST_VECT_HPP__

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

//modif AL AL
// modif AL pour ajouter la fonction gemm
#include "faust_linear_algebra.h"

template<typename FPP>
const char * Faust::Vect<FPP,Cpu>::class_name = "Faust::Vect<FPP,Cpu>::";

template<typename FPP>
template<typename FPP1>
Faust::Vect<FPP,Cpu>::Vect(const Faust::Vect<FPP1,Cpu>& v) : dim(v.dim), vec(v.dim)
{
for (int i=0;i<dim;i++)
	vec[i]=v(i);
}


template<typename FPP>
bool Faust::Vect<FPP,Cpu>::equality(Faust::Vect<FPP,Cpu> const &x, FPP precision) const
{
  if (size() != x.size())
	handleError(class_name," equality : different dimension");

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
Faust::Vect<FPP,Cpu>::Vect(const faust_unsigned_int dim_, const FPP* data_) : dim(dim_), vec(dim_)
{
		memcpy(getData(), data_, dim*sizeof(FPP));
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
       cout << "empty vector";
    for (int i=0 ; i<size() ; i++)
	    cout << getData()[i] << " " ;
    cout<<endl<<endl;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::print_file(const char* filename)const
{
	ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<size() ; i++)
		fichier << setprecision(20) << vec (i) << endl;
	fichier.close();
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::resize(const int new_dim)
{

		if (new_dim <0)
		{
			handleError(class_name,"resize : new dimensions must be positive");
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
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha;
}

template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator+=(const Faust::Vect<FPP,Cpu>& v)
{
   if(v.size()!=size())
   {

	  handleError(class_name,"operator+= : dimensions are in conflict");
   }
   FPP*const ptr_data = getData();
   FPP*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += v_ptr_data[i];
}


template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator-=(const Faust::Vect<FPP,Cpu>& v)
{
   if(v.size()!=size())
   {
	   handleError(class_name,"operator-= : dimensions are in conflict");
   }
   FPP*const ptr_data = getData();
   FPP*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] -= v_ptr_data[i];
}



template<typename FPP>
FPP Faust::Vect<FPP,Cpu>::mean_relative_error(const Faust::Vect<FPP,Cpu>& v_ref)
{
   if(v_ref.size() != size())
   {
     handleError(class_name,"relative_error : sizes are different");

   }

   Faust::Vect<FPP,Cpu> tmp(size());
   for(int i=0 ; i<size() ; i++)
      tmp[i] = fabs((vec[i]-v_ref.vec[i])/v_ref.vec[i]);

   return tmp.mean();
}


template<typename FPP>
template<typename FPP1>
void Faust::Vect<FPP,Cpu>::operator=(Faust::Vect<FPP1,Cpu> const& y)
{

	resize(y.dim);
	for (int i=0;i<dim;i++)
		vec[i]=y(i);

}


template<typename FPP>
void Faust::Vect<FPP,Cpu>::operator=(Faust::Vect<FPP,Cpu> const& y)
{
	  vec = y.vec;
	  dim = y.dim;
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
void  Faust::Vect<FPP,Cpu>::multiplyLeft(Faust::MatSparse<FPP,Cpu> const& A)
{
	int nbColA = A.getNbCol();

	if(nbColA != vec.size())
	{
		 handleError(class_name,"multiplyLeft : incorrect dimensions");

	}

		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
		#endif
		vec = A.mat * vec;
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			cout <<"0 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			t_local_multiplyLeft.reset();
		#endif


	dim = A.dim1;
}

#endif
