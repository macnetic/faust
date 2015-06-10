#include "faust_init_from_matio.h"
#include <iostream>
#include "faust_constant.h"
#include "faust_mat.h"

using namespace std;


int main()
{
   faust_mat A(10,10);
   faust_real B,E,G;
   int C;
   bool D,F;

   init_faust_mat_from_matio_mat(A, "test.mat", "A");
   B = init_faust_mat_from_matio_double("test.mat", "B");
   C = init_faust_mat_from_matio_int("test.mat", "C");
   D = init_faust_mat_from_matio_bool("test.mat", "D");
   F = init_faust_mat_from_matio_bool("test.mat", "F");
   //G = init_faust_mat_from_matio_double("test.mat", "G");


  cout<<"A="<<endl;
  for (int i=0 ; i<A.getNbRow() ; i++)
  { 
    for (int j=0 ; j<A.getNbCol() ;j++)
     {
        cout << A.getData()[j*A.getNbRow()+i]<<" ";
     }
     cout<<endl;
  }

  cout << endl << "B="<<B<<endl;
  cout << endl << "C="<<C<<endl;
  cout << endl << "D="<<D<<endl;
  cout << endl << "F="<<F<<endl;
   
   


return 0;
}


