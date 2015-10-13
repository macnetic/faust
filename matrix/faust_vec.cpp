#include "faust_vec.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "faust_mat.h"
#include "faust_spmat.h"
#include "faust_exception.h"

using namespace std;


#ifdef __GEMM_WITH_MKL__
	#include "mkl_spblas.h"
#endif

const char * faust_vec::class_name = "faust_vec::";

faust_vec::faust_vec(const faust_unsigned_int dim_, const faust_real* data_) : dim(dim_), vec(dim_)
{
		memcpy(getData(), data_, dim*sizeof(faust_real));	
}

void faust_vec::setOnes()
{
	vec.setOnes();
}


void faust_vec::Display() const
{
	 cout << "dim = " << size() << endl;
	 cout << vec << endl;
}
 
void faust_vec::print_file(const char* filename)const
{
	ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<size() ; i++)
		fichier << setprecision(20) << vec (i) << endl;
	fichier.close();
}


void faust_vec::resize(const int new_dim)
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

void faust_vec::operator*=(const faust_real alpha)
{
   faust_real*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] *= alpha; 
}
void faust_vec::operator+=(const faust_real alpha)
{
   faust_real*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha; 
}
void faust_vec::operator-=(const faust_real alpha)
{
   faust_real*const ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += alpha; 
}


void faust_vec::operator+=(const faust_vec& v)
{
   if(v.size()!=size())
   {

	  handleError(class_name,"operator+= : dimensions are in conflict");
   }
   faust_real*const ptr_data = getData();
   faust_real*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] += v_ptr_data[i]; 
}

void faust_vec::operator-=(const faust_vec& v)
{
   if(v.size()!=size())
   {
	   handleError(class_name,"operator-= : dimensions are in conflict");
   }
   faust_real*const ptr_data = getData();
   faust_real*const v_ptr_data = getData();
   for (int i=0 ; i<size() ; i++)
      ptr_data[i] -= v_ptr_data[i]; 
}

faust_real faust_vec::mean_relative_error(const faust_vec& v_ref)
{
   if(v_ref.size() != size())
   {
     handleError(class_name,"relative_error : sizes are different");

   }

   faust_vec tmp(size());
   for(int i=0 ; i<size() ; i++)
      tmp[i] = fabs((vec[i]-v_ref.vec[i])/v_ref.vec[i]);

   return tmp.mean();
}




void faust_vec::operator=(faust_vec const& y)
{
	  vec = y.vec;
	  dim = y.dim;

}


void  faust_vec::multiplyLeft(faust_spmat const& A)
{	
	int nbColA = A.getNbCol();
	
	if(nbColA != vec.size())
	{
		 handleError(class_name,"multiplyLeft : incorrect dimensions");
		
	}
	#ifdef __GEMM_WITH_MKL__
		int nbRowA = A.getNbRow();
		faust_real alpha = 1.0,beta = 0.0;
		char matdescrA[8];
		matdescrA[0]='G';
		for (int i=1 ; i<=6 ;i++)
		matdescrA[i] = 'C';
		matdescrA[7] = '\0';
		faust_vec y(nbRowA);
		faust_spmat Acopy(A);
		if (!Acopy.isCompressedMode())
		{
			Acopy.makeCompression();
			std::cout<<class_name<<"multiplyLeft : sparse matrix A is not compressed (possible lack of time)"<<std::endl;
		}
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
		#endif	
		#ifdef FAUST_SINGLE	
			mkl_scsrmv("T",&nbColA,&nbRowA,&alpha,matdescrA,Acopy.getValuePtr(),Acopy.getInnerIndexPtr(),Acopy.getOuterIndexPtr(),&(Acopy.getOuterIndexPtr()[1]),getData(),&beta,y.getData());
		#else
			//cout << "MKL debut" << endl;
			mkl_dcsrmv("T",&nbColA,&nbRowA,&alpha,matdescrA,Acopy.getValuePtr(),Acopy.getInnerIndexPtr(),Acopy.getOuterIndexPtr(),&(Acopy.getOuterIndexPtr()[1]),getData(),&beta,y.getData());
			//cout << "MKL fin" << endl;
		#endif
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			cout << "1 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			t_local_multiplyLeft.reset();
		#endif	
		(*this) = y;	
	#else
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
		#endif	
		vec = A.mat * vec;
		#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			cout <<"0 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			t_local_multiplyLeft.reset();
		#endif		
	#endif
	
	dim = A.dim1;
}


