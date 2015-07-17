#include "faust_vec.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "faust_mat.h"
#include "faust_spmat.h"

using namespace std;

void faust_vec::setOnes()
{
	vec.setOnes();
}


void faust_vec::Display() const
{
	 cout << "dim = " << getDim() << endl;
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
			cerr << "ERROR FAUST_VEC resize : les nouvelles dimensions doivent etre strictement positive" << endl;
			exit( EXIT_FAILURE);
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
      cerr << "Error in faust_vec::operator+= : sizes are different" << endl;
      exit(EXIT_FAILURE);
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
      cerr << "Error in faust_vec::operator-= : sizes are different" << endl;
      exit(EXIT_FAILURE);
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
      cerr << "Error in faust_vec::relative_error : sizes are different" << endl;
      exit(EXIT_FAILURE);
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
		cerr << "Error in faust_vec::multiplyLeft : incorrect dimensions" << endl;
		exit(EXIT_FAILURE);
	}
	vec = A.mat * vec;
	dim = A.dim1;
}


