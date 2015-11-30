#include <iostream>
#include <iomanip>
#include <fstream>


using namespace std;


#ifdef __GEMM_WITH_MKL__
	#include "mkl_spblas.h"
#endif
template<typename T>
const char * faust_vec<T>::class_name = "faust_vec<T>::";

template<typename T>
template<typename U>
faust_vec<T>::faust_vec(const faust_vec<U>& v) : dim(v.dim), vec(v.dim)
{
for (int i=0;i<dim;i++)
	vec[i]=v(i);	
}

template<typename T>
faust_vec<T>::faust_vec(const faust_unsigned_int dim_, const T* data_) : dim(dim_), vec(dim_)
{
		memcpy(getData(), data_, dim*sizeof(T));	
}

template<typename T>
void faust_vec<T>::setOnes()
{
	vec.setOnes();
}

template<typename T>
void faust_vec<T>::Display() const
{
	 cout << "dim = " << size() << endl;
	 cout << vec << endl;
}
 
template<typename T> 
void faust_vec<T>::print_file(const char* filename)const
{
	ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<size() ; i++)
		fichier << setprecision(20) << vec (i) << endl;
	fichier.close();
}

template<typename T>
void faust_vec<T>::resize(const int new_dim)
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

template<typename T>
void faust_vec<T>::operator*=(const T alpha)
{
   T*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] *= alpha; 
}

template<typename T>
void faust_vec<T>::operator+=(const T alpha)
{
   T*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha; 
}


template<typename T>
void faust_vec<T>::operator-=(const T alpha)
{
   T*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha; 
}

template<typename T>
void faust_vec<T>::operator+=(const faust_vec<T>& v)
{
   if(v.size()!=size())
   {

	  handleError(class_name,"operator+= : dimensions are in conflict");
   }
   T*const ptr_data = getData();
   T*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += v_ptr_data[i]; 
}


template<typename T>
void faust_vec<T>::operator-=(const faust_vec<T>& v)
{
   if(v.size()!=size())
   {
	   handleError(class_name,"operator-= : dimensions are in conflict");
   }
   T*const ptr_data = getData();
   T*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] -= v_ptr_data[i]; 
}



template<typename T>
T faust_vec<T>::mean_relative_error(const faust_vec<T>& v_ref)
{
   if(v_ref.size() != size())
   {
     handleError(class_name,"relative_error : sizes are different");

   }

   faust_vec<T> tmp(size());
   for(int i=0 ; i<size() ; i++)
      tmp[i] = fabs((vec[i]-v_ref.vec[i])/v_ref.vec[i]);

   return tmp.mean();
}


template<typename T>
template<typename U>
void faust_vec<T>::operator=(faust_vec<U> const& y)
{	
	  
	resize(y.dim);
	for (int i=0;i<dim;i++)
		vec[i]=y(i);

}


template<typename T>
void faust_vec<T>::operator=(faust_vec<T> const& y)
{
	  vec = y.vec;
	  dim = y.dim;

}

template<typename T>
void  faust_vec<T>::multiplyLeft(faust_spmat<T> const& A)
{	
	int nbColA = A.getNbCol();
	
	if(nbColA != vec.size())
	{
		 handleError(class_name,"multiplyLeft : incorrect dimensions");
		
	}
	// #ifdef __GEMM_WITH_MKL__
		// int nbRowA = A.getNbRow();
		// T alpha = 1.0,beta = 0.0;
		// char matdescrA[8];
		// matdescrA[0]='G';
		// for (int i=1 ; i<=6 ;i++)
		// matdescrA[i] = 'C';
		// matdescrA[7] = '\0';
		// faust_vec<T> y(nbRowA);
		// faust_spmat<T> Acopy(A);
		// if (!Acopy.isCompressedMode())
		// {
			// Acopy.makeCompression();
			// std::cout<<class_name<<"multiplyLeft : sparse matrix A is not compressed (possible lack of time)"<<std::endl;
		// }
		// #ifdef __COMPILE_TIMERS__
			// t_local_multiplyLeft.start();
		// #endif	
		// #ifdef FAUST_SINGLE	
			// mkl_scsrmv("T",&nbColA,&nbRowA,&alpha,matdescrA,Acopy.getValuePtr(),Acopy.getInnerIndexPtr(),Acopy.getOuterIndexPtr(),&(Acopy.getOuterIndexPtr()[1]),getData(),&beta,y.getData());
		// #else
			// cout << "MKL debut" << endl;
			// mkl_dcsrmv("T",&nbColA,&nbRowA,&alpha,matdescrA,Acopy.getValuePtr(),Acopy.getInnerIndexPtr(),Acopy.getOuterIndexPtr(),&(Acopy.getOuterIndexPtr()[1]),getData(),&beta,y.getData());
			// cout << "MKL fin" << endl;
		// #endif
		// #ifdef __COMPILE_TIMERS__
			// t_local_multiplyLeft.stop();
			// cout << "1 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			// t_local_multiplyLeft.reset();
		// #endif	
		// (*this) = y;	
	// #else
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
		#endif	
		vec = A.mat * vec;
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			cout <<"0 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			t_local_multiplyLeft.reset();
		#endif		
	// #endif
	
	dim = A.dim1;
}


