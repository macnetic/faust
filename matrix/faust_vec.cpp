#include "faust_vec.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "faust_mat.h"
#include "faust_spmat.h"
#include "faust_exception.h"

using namespace std;

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
	if(A.getNbCol() != vec.size())
	{
		 handleError(class_name,"multiplyLeft : incorrect dimensions");
		
	}
	vec = A.mat * vec;
	dim = A.dim1;
}


